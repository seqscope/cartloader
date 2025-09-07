```bash
zenodo_token=/path/to/zenodo/token/file    # replace /path/to/zenodo/token/file by path to your zenodo token file

docker run -it --rm \
  -v $(pwd):/data \
  weiqiuc/cartloader:20250708b \
  upload_zenodo \
    --in-dir /data/cartload2 \
    --upload-method catalog \
    --zenodo-token ${zenodo_token} \
    --create-new-deposition \
    --title  "Yur Title" \
    --creators "Your Name" \
    --description "This is an example description"
```