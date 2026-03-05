```bash
zenodo_token=/path/to/zenodo/token/file    # replace /path/to/zenodo/token/file with the path to your Zenodo token file

docker run -it --rm \
  -v $(pwd):/data \
  weiqiuc/cartloader:${docker_tag} \
  upload_zenodo \
    --in-dir /data/cartload2 \
    --upload-method catalog \
    --zenodo-token ${zenodo_token} \
    --title  "Your Title" \
    --creators "Your Name" \
    --description "This is an example description"
```
