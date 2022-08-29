## Issues

- Unexpected movements are generated.
  - ~~`unif_index2()` correct~~
  - `prof_cells` unexpected
    - number of cells in `prof_cells`: match
    - cells values: match but in wrong order
      - order of `prof_cells` is controlled by `cell_den`
      - order of `cell_den` is expected
      - `unif_index` generates samples in wrong order

  Solution: Stop struggling with reproducing extactly same R behaviour.