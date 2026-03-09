# What are the Zenodo token and deposition ID, and how do you get them?

## Zenodo Token File

`CartLoader` uses the Zenodo API for uploads. To authenticate, Zenodo requires an access token.

To obtain a token:

1. Log in to [Zenodo](https://zenodo.org).
2. Open [Applications](https://zenodo.org/account/settings/applications/).
3. Click "New Token" and select scopes such as `deposit:write` and `deposit:actions`.
4. Copy the generated token.
5. Save it in a plain text file and pass that file path to `--zenodo-token`.

## Zenodo Deposition ID

A deposition ID is the numeric identifier of a Zenodo deposition (dataset record).  
If you already created one, find it in the URL, for example:

```text
https://zenodo.org/deposit/1234567
                            ^^^^^^^
                           deposition ID
```

See also:

- [Zenodo upload reference](../reference/upload_zenodo.md)
