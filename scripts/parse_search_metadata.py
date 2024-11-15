import os
import pandas as pd

MONTHS = {
    "JAN": "01",
    "FEB": "02",
    "MAR": "03",
    "APR": "04",
    "MAY": "05",
    "JUN": "06",
    "JUL": "07",
    "AUG": "08",
    "SEP": "09",
    "OCT": "10",
    "NOV": "11",
    "DEC": "12",
}


def parse_metadata(sample):

    if "ENC" in sample:
        loc = "Encina"
        date = sample[sample.index(abbrev) + 3 : sample.index("ENC") + 8]
    elif "SB" in sample:
        abbrev = "SB"
        loc = "South Bay"
        date = sample[sample.index(abbrev) + 2 : sample.index("SB") + 7]
    elif "PL" in sample:
        abbrev = "PL"
        loc = "Point Loma"
        date = sample[sample.index(abbrev) + 2 : sample.index("PL") + 7]

    month = MONTHS[date[:3]]
    day = date[3:]

    submission_date = sample[: sample.index(abbrev)]
    submission_month = submission_date[:2]
    submission_year = submission_date[-3:-1]

    if submission_month == "01" and month == "DEC":
        year = f"20{int(submission_year)-1}"
    else:
        year = f"20{submission_year}"

    collection_date = f"{year}-{month}-{day}"

    return collection_date, loc

def main():
    metadata = pd.DataFrame(columns=["Sample", "Collection_Date", "Location"])

    metadata["Sample"] = [f.split(".trimmed")[0] for f in os.listdir("bam")]

    metadata["Collection_Date"], metadata["Location"] = zip(
        *metadata["Sample"].map(parse_metadata)
    )
    metadata["Collection_Date"] = pd.to_datetime(metadata["Collection_Date"])
    metadata = metadata.sort_values("Collection_Date")
    metadata.to_csv("SEARCH_metadata.csv", index=False)

if __name__ == '__main__':
    main()