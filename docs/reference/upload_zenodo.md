# Upload to Zenodo

!!! tip "Not part of `run_together`"
    [`run_together`](./run_together.md#publishing) has no Zenodo stage; its publishing step uploads to S3 only. Run `upload_zenodo` yourself on a packaged directory (e.g., `cartl/<id>/`).

## Overview

Use `upload_zenodo` to publish `CartLoader` outputs (PMTiles, factor files, and the catalog) to a Zenodo deposition, either an existing one or a newly created draft. In catalog mode, files are taken from `catalog.yaml` and uploaded in three groups: cartload basics (with `catalog.yaml`), additional basemaps (non-SGE PMTiles), and optional files (UMAP, alias), using the same grouping as [`upload_aws`](./upload_aws.md).

!!! danger "`--restart` deletes every file in the target draft"
    With `--zenodo-deposition-id`, `--restart` removes **all** files currently in the target draft deposition before uploading, including files that are not part of this upload and files a new-version draft inherited from the published record. To replace only the files you are uploading, use `--overwrite` instead.

!!! warning "`--dry-run` still changes Zenodo"
    `--dry-run` skips file uploads, file deletions, and publishing, but it still:

    - creates a new deposition when `--zenodo-deposition-id` is omitted;
    - opens a new-version draft when the ID points to a published record;
    - writes the merged metadata to the deposition.

---
## Requirements

- A packaged directory from [`run_cartload2`](./run_cartload2.md) with a `catalog.yaml` (for `--upload-method catalog`).
- A Zenodo access token saved in a file.
- For a new deposition: `--title` and at least one `--creators` entry.
- (Optional) A Zenodo deposition ID to upload to an existing deposition.

!!! question
    - [What are the Zenodo token and deposition ID, and how do you get them?](../faq/zenodo_setup.md)

---
## Actions

1. **Choose the deposition.** Without `--zenodo-deposition-id`, a new draft is created with the given metadata. With an ID of a draft, that draft is used. With an ID of a published record, its latest unsubmitted draft is reused if one exists; otherwise a new version is created. New-version drafts start with the files of the published version.
2. **Update metadata** (existing depositions). Provided `--title`, `--description`, and `--creators` replace the existing values; `--upload-type` (default `dataset`) always replaces the existing type.
3. **Select files** by `--upload-method`. Files are uploaded under their basenames.
4. **Check capacity.** If the deposition would hold more than 100 files (existing plus new distinct names), the command aborts before uploading.
5. **Upload.** Files not yet in the deposition are uploaded; files whose names already exist are skipped unless `--overwrite` (or `--restart`) is set. In catalog mode the groups are uploaded in order: basics with `catalog.yaml`, then basemaps, then optional files.
6. **Publish** with `--publish` (skipped in dry-run mode). Uploads rejected by Zenodo are listed at the end with a command to retry them; other upload errors (e.g., connection failures) are only printed in the log.

---
## Example Usage

=== "Create a new deposition to upload"

    ```bash
    cartloader upload_zenodo \
      --upload-method catalog \
      --in-dir /path/to/run_cartload2/output/directory \
      --zenodo-token /path/to/zenodo_token.txt \
      --title  "Title Info" \
      --creators "Lastname, Firstname" \
      --description "Description Info"
    ```

    - `--upload-method catalog`: Upload the files referenced by the catalog.
    - `--title` and `--creators` are required for a new deposition.

=== "Upload to an existing deposition"

    ```bash
    zenodo_deposition_id=DEPOSITION_ID                # replace with your deposition ID (published IDs get a new draft automatically)

    cartloader upload_zenodo \
      --upload-method catalog \
      --in-dir /path/to/run_cartload2/output/directory \
      --zenodo-token /path/to/zenodo_token.txt \
      --zenodo-deposition-id ${zenodo_deposition_id}
    ```

    Add `--overwrite` to replace files that already exist in the deposition.

=== "Upload explicit files"

    ```bash
    cartloader upload_zenodo \
      --upload-method files \
      --in-dir /path/to/run_cartload2/output/directory \
      --files catalog.yaml genes_all.pmtiles some_factor.pmtiles \
      --zenodo-token /path/to/zenodo_token.txt \
      --zenodo-deposition-id ${zenodo_deposition_id}
    ```

    - `--upload-method files`: Upload only the files given by `--files`.

---

## Parameters

### Input Selection
- `--in-dir` (str, required): Directory containing the files to upload.
- `--upload-method` (str): Which files to upload (default: `all`).

    !!! info "`--upload-method` Options"
        * `all`: every regular file directly in `--in-dir` (not recursive; includes Makefiles, logs, and flag files);
        * `catalog`: the files referenced by `catalog.yaml`, plus the catalog itself;
        * `files`: only the files given by `--files`.

- `--files` (list of str): Files to upload, relative to `--in-dir` or absolute. Required with `--upload-method files`.
- `--catalog-yaml` (str): Catalog for `--upload-method catalog` (default: `<in-dir>/catalog.yaml`).

### Zenodo Configuration
- `--zenodo-token` (str, required): Path to a file containing your Zenodo access token.
- `--zenodo-deposition-id` (str): Existing deposition ID. A published record gets a draft (reused or newly created); omit to create a new deposition.
- `--publish` (flag): Publish the deposition after upload. Recommended to leave off and publish manually after reviewing the draft on Zenodo.

### Zenodo Metadata
- `--title` (str): Deposition title. Required for a new deposition.
- `--creators` (list of str): Creators, each quoted, e.g. `"Lastname, Firstname"`. Required for a new deposition.
- `--description` (str): Deposition description (new depositions default to `No description provided.`).
- `--upload-type` (str): One of `dataset`, `software`, `publication`, `poster`, `presentation`, `image`, `video`, `lesson`, `other` (default: `dataset`; also applied to existing depositions).

### Behavior Flags
- `--overwrite` (flag): Re-upload files whose names already exist in the deposition, replacing them.
- `--restart` (flag): Delete **all** files in the target draft first, then upload (see the warning above). Implies `--overwrite`.
- `--dry-run` (flag): Print the files that would be uploaded or deleted; Zenodo is still contacted and modified as described above.

---
## Output

- Files uploaded to the Zenodo deposition (existing draft, new version draft, or new deposition); the deposition ID is printed.
- With `--upload-method catalog` (not in dry-run mode): `<in-dir>/cartload.zenodo.done` is touched after the basics group, even if some of its files failed to upload, and `<file>.zenodo.done` is written next to each successfully uploaded basemap or optional file. Check the log for failures before relying on these flags.
