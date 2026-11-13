# Package configuration

`stratonix package <description.yaml> -o <new folder> --config package.yaml [--zip]`

All keys are optional.

```yaml
formats: [petrel, vtk, csv]       # export ids or sets; default [petrel]
figures: true                     # figure set in SVG, PNG, TIFF and PPTX; default true
zip: true                         # also write <folder>.zip; default false (or --zip)
terms_of_use: "..."               # text approved by the model owner; without it the package is a DRAFT
contact: "..."                    # overrides model.contact of the description
assumptions: ["..."]              # added to the README (assumed parameters are listed automatically)
limitations: ["..."]              # added to the README
banned_patterns: ['PRJ-[0-9]{4}'] # regular expressions that must not appear in any text file or file name
translations:                     # internal label -> external wording, applied to every text file
  "PRJ-0001": "the client project"
default_patterns: true            # also apply the built-in patterns below
```

Built-in banned patterns (Stratonix defaults): register identifiers (`D-001`, `I-004`, `O-010`,
`S-002`, `R-1.0.0`), work-breakdown codes (`WBS-...`), internal note identifiers (`STX-TN-..`,
`STX-RI-..`) and absolute local paths (`C:\...`, `/home/...`, `/Users/...`, `/mnt/...`, `/tmp/...`).

Order of work: inputs and their sharing status are listed and an input with status `unknown`
stops the command; the model is built and checked (a failed check stops it); files are written;
translations are applied; the scan runs. If the scan finds anything, the folder is renamed
`<folder>_REJECTED_SCAN` with `SCAN_FAILED.txt` listing every hit, and no zip is written. Otherwise
`manifest.json` (SHA-256 of every file) and the optional zip are written. The manifest is
portable: it records the Stratonix version and commit but not the compiler or platform, so the
same package built on two machines differs only by the manifest timestamp (raster figures can
also differ where the machines have different fonts).
