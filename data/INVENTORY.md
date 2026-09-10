# Nature-style source-data workbooks

These builders live in `data/` with the input tables. They are optional and are
not required to regenerate the manuscript figures.

`build_source_data.py` writes one workbook per figure into `data/` **after** the
analysis scripts have been run. Those Excel files are generated output and are
not committed.

```bash
python data/build_source_data.py
```

`build_negative_scores_supplementary_table.py` summarises negative DDA scores
across the full compound–profile matrix. That matrix is not redistributed with
this repository. The script therefore requires an explicit path:

```bash
python data/build_negative_scores_supplementary_table.py --compounds-csv PATH/TO/score_matrix.csv
```
