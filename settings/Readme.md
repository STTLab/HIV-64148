# Settings

This directory contains the following files.

1. settings.json
2. secrets.json

The files are in JSON (key:value) format.
You can modify values within these files according to your preference.
Please note that alteration of key name or value's datatype may resulted in programmatic error.

## secrets.json

### API keys

#### NCBI - optional (str)

```JSON
"API_KEYS": {
    ...
    "NCBI": "____YOUR_KEY____"
    ...
}
```

Use to increase request rate to E-utils API. Users are allowed 3 requests/second without an API key. Create an API key to increase your e-utils limit to 10 requests/second. The key can be obtained by singing up at [NCBI website](https://www.ncbi.nlm.nih.go) via the [account settings](https://www.ncbi.nlm.nih.gov/account/settings/) page.
