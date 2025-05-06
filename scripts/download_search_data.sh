set -ex

while read -r line; do
  aws s3 cp s3://ucsd-all/${line} san_diego --recursive --exclude "*" --include "*.trimmed.sorted*.bam*"
  aws s3 cp s3://ucsd-all/${line} raw_metadata --recursive --exclude "*" --include "*.csv"

done < data/search_bucket_2023-2025.txt
