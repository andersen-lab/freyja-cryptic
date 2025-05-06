import pandas as pd
import os 


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
        date = sample[sample.index("ENC") + 3 : sample.index("ENC") + 8]
        submission_date = sample[: sample.index("ENC")]
    elif "SB" in sample:
        loc = "South Bay"
        date = sample[sample.index("SB") + 2 : sample.index("SB") + 7]
        submission_date = sample[: sample.index("SB")]
    elif "PL" in sample:
        loc = "Point Loma"
        date = sample[sample.index("PL") + 2 : sample.index("PL") + 7]
        submission_date = sample[: sample.index("PL")]

    try:
        month = MONTHS[date[:3]]
        day = date[3:]
    except KeyError:
        print(f"Error: {sample}")
        return None, None

    
    submission_month = submission_date[:2]
    submission_year = submission_date[-3:-1]

    if submission_month == "01" and month == "12":
        year = f"20{int(submission_year)-1}"
    else:
        year = f"20{submission_year}"

    collection_date = f"{year}-{month}-{day}"

    return collection_date, loc

def main():
    META_DIR = "metadata/"

    combined_metadata = pd.DataFrame()
    
    # Read all metadata files
    for file in os.listdir(META_DIR):
        if file.endswith(".csv"):
            file_path = os.path.join(META_DIR, file)
            metadata = pd.read_csv(file_path)
            
            # Combine metadata
            combined_metadata = pd.concat([combined_metadata, metadata], ignore_index=True)

    print("Combined_metadata.columns:", combined_metadata.columns)
    combined_metadata["collection_date"], combined_metadata["location"] = zip(
        *combined_metadata["Sample"].map(parse_metadata)
    )
    combined_metadata["collection_date"] = pd.to_datetime(combined_metadata["collection_date"])
    combined_metadata = combined_metadata.sort_values("collection_date")

    # Save combined metadata
    combined_metadata.to_csv("combined_metadata.csv", index=False)

if __name__ == '__main__':
    main()