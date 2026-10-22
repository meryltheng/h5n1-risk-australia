"""
Download the "Events by species" table behind the DAFF bird flu page.

https://www.agriculture.gov.au/campaigns/birdflu/latest-data

The table is a Power BI "publish to web" embed, so there is no HTML table to
scrape. Instead this queries the same public Power BI API the embed uses, then
decodes Power BI's compressed response into a table.

The result is a snapshot of the table on the day the script is run. It cannot
retrieve the table as it was on an earlier date.
"""

import datetime
import json

import pandas as pd
import requests

# Resource key of the "Events by species" embed: the "k" value inside the
# base64 `r=` parameter of the iframe under the "Events by species" heading.
# If DAFF republishes the report, this key changes.
RESOURCE_KEY = "efa7acc9-75fa-48a9-93f6-6e34c9348df2"
API_BASE = "https://wabi-australia-southeast-api.analysis.windows.net/public/reports"
HEADERS = {"X-PowerBI-ResourceKey": RESOURCE_KEY}
MAX_ROWS = 30_000  # upper limit per request; the table currently has ~600

COLUMNS = [
    "date_sampled",
    "jurisdiction",
    "lga",
    "locality",
    "common_name",
    "scientific_name",
    "record_id",
]


def col_ref(source: str, prop: str) -> dict:
    """Power BI's way of referring to a column of a table."""
    return {"Column": {"Expression": {"SourceRef": {"Source": source}}, "Property": prop}}


def build_query_body(model_id, dataset_id, report_id, visual_id) -> dict:
    """The same columns, sort order and filters as the table on the DAFF page:
    the latest complete reporting date, and jurisdiction not 'EXT' (external
    territories) or blank."""
    query = {
        "Version": 2,
        "From": [
            {"Name": "p", "Entity": "Positive events", "Type": 0},
            {"Name": "r", "Entity": "Reporting calendar", "Type": 0},
            {"Name": "j", "Entity": "_Jurisdiction", "Type": 0},
        ],
        "Select": [
            {**col_ref("p", col), "Name": f"Positive events.{col}"} for col in COLUMNS
        ],
        "Where": [
            {"Condition": {"In": {
                "Expressions": [col_ref("r", "Is latest complete date (Event summary)")],
                "Values": [[{"Literal": {"Value": "true"}}]],
            }}},
            {"Condition": {"Not": {"Expression": {"In": {
                "Expressions": [col_ref("j", "Jurisdiction")],
                "Values": [
                    [{"Literal": {"Value": "'EXT'"}}],
                    [{"Literal": {"Value": "null"}}],
                ],
            }}}}},
        ],
        "OrderBy": [{"Direction": 2, "Expression": col_ref("p", "date_sampled")}],
    }

    return {
        "version": "1.0.0",
        "queries": [{
            "Query": {"Commands": [{"SemanticQueryDataShapeCommand": {
                "Query": query,
                "Binding": {
                    "Primary": {"Groupings": [{"Projections": list(range(len(COLUMNS)))}]},
                    "DataReduction": {
                        "DataVolume": 3,
                        "Primary": {"Window": {"Count": MAX_ROWS}},
                    },
                    "Version": 1,
                },
            }}]},
            "QueryId": "",
            "ApplicationContext": {
                "DatasetId": dataset_id,
                "Sources": [{"ReportId": report_id, "VisualId": visual_id}],
            },
        }],
        "cancelQueries": [],
        "modelId": model_id,
    }


def decode_rows(rows: list[dict], dicts: dict) -> list[list]:
    """Decode Power BI's compressed rows into plain lists of values.

    Each row lists only the values that changed, in "C":
      R = bitmask of columns that repeat the previous row's value
      Ø = bitmask of columns that are null
    Text values are usually indices into a lookup ("ValueDicts"); once a lookup
    reaches its cap, further values arrive as literal strings instead.
    The column schema ("S") is given on the first row only.
    """
    schema = rows[0]["S"]
    previous = [None] * len(schema)
    decoded = []

    for row in rows:
        repeat_mask = row.get("R", 0)
        null_mask = row.get("Ø", 0)
        new_values = iter(row.get("C", []))
        current = []

        for i, column in enumerate(schema):
            bit = 1 << i  # the bit for column i
            if repeat_mask & bit:
                value = previous[i]
            elif null_mask & bit:
                value = None
            else:
                value = next(new_values)
                dict_name = column.get("DN")
                if dict_name is not None and isinstance(value, int):
                    value = dicts[dict_name][value]  # indices are 0-based, as in Python
            current.append(value)

        decoded.append(current)
        previous = current

    return decoded


def run(output_path: str, snapshot_date: str) -> None:
    today = datetime.date.today()
    if snapshot_date != today.strftime("%Y%m%d"):
        raise ValueError(
            f"Requested a DAFF snapshot for {snapshot_date}, but this script can "
            "only download today's table."
        )

    # 1. Report metadata: dataset, report and visual IDs
    response = requests.get(
        f"{API_BASE}/{RESOURCE_KEY}/modelsAndExploration",
        params={"preferReadOnlySession": "true"},
        headers=HEADERS,
        timeout=60,
    )
    response.raise_for_status()
    meta = response.json()

    model = meta["models"][0]
    report_id = meta["exploration"]["report"]["objectId"]

    # find the table visual ("tableEx") in the report definition
    doc = json.loads(meta["exploration"]["explorationContent"]["explorationDocument"])
    visual_id = None
    for page in doc["pages"]["pages"]:
        for container in page.get("visualContainers", []):
            content = container.get("content", {})
            if content.get("visual", {}).get("visualType") == "tableEx":
                visual_id = content["name"]
    if visual_id is None:
        raise RuntimeError("No table visual found; the report layout has changed.")

    # 2. Query the table
    body = build_query_body(model["id"], model["dbName"], report_id, visual_id)
    response = requests.post(
        f"{API_BASE}/querydata",
        params={"synchronous": "true"},
        json=body,
        headers=HEADERS,
        timeout=60,
    )
    response.raise_for_status()
    ds = response.json()["results"][0]["result"]["data"]["dsr"]["DS"][0]

    # 3. Decode into a table
    rows = ds["PH"][0]["DM0"]
    events = pd.DataFrame(decode_rows(rows, ds.get("ValueDicts", {})), columns=COLUMNS)

    # date_sampled arrives as milliseconds since 1 January 1970 (UTC)
    events["date_sampled"] = pd.to_datetime(
        pd.to_numeric(events["date_sampled"]), unit="ms", utc=True
    ).dt.strftime("%Y-%m-%d")
    events["data_refreshed_utc"] = model["LastRefreshTime"]
    events["retrieved_on"] = today.isoformat()

    # 4. Stop rather than save an incomplete table
    if ds.get("IC") is not True or "RT" in ds or len(events) >= MAX_ROWS:
        raise RuntimeError("Result may be truncated; increase MAX_ROWS.")

    events.to_csv(output_path, index=False)
    print(f"{len(events)} events; {events['date_sampled'].isna().sum()} with no date sampled.")


if __name__ == "__main__":
    run(
        output_path=snakemake.output.csv, 
        snapshot_date=snakemake.wildcards.date
        )