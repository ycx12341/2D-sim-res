## Issues

- ~~Unexpected movements are generated.~~
  - ~~`unif_index2()` correct~~
  - ~~`prof_cells` unexpected~~
    - ~~number of cells in `prof_cells`: match~~
    - ~~cells values: match but in wrong order~~
      - ~~order of `prof_cells` is controlled by `cell_den`~~
      - ~~order of `cell_den` is expected~~
      - ~~`unif_index` generates samples in wrong order~~

  Solution: Stop struggling with reproducing extactly same R behaviour.

- ~~Can not properly fetch data in range 
  <br> `ind_pos[y.cut[ii]:(y.cut[ii] + (mat.size - 1)), x.cut[jj]:(x.cut[jj] + (mat.size - 1))]`~~
  - ~~`ind_pos`: correct~~
  - ~~`iter_range()`: incorrect~~
    - ~~using `i` instead of `ii`~~
    - ~~Solved~~

- Program crashed in `generate_pattern()`
  - seed: 234567
  - dimension: 243-th
  - on call: `initial_condition()`
  - `INIT_CELLS_COLS` in `Parameters`: error
    - total dimensions: 500
      ```
      3.92764 2.26746 2.23225 4.3747 1.93549 4.53162 1.79133 4.56243 1.54211 2.7872 0 1.29741e-231 0.0111366 0.0926829 0.0673112 
      0.0394321 0.0797405 0.0447875 0.0520084 0.0577969 0.075692 0.0621816 0 1.29741e-231 0.476869 0.498841 0.437558 
      0.391151 0.364583 0.52035 0.257124 0.769458 0.24769 0.222052 0 1.29741e-231 3.76968e-317 0 0 0 0 0 0 0 0 0 0 0 0 
      6.95259e-310 0 0 0 6.95259e-310 0 0 0 6.95259e-310 0 0 0 6.95259e-310 0 0 0 6.95259e-310 0 0 0 2.122e-314 0 0 0 
      6.95259e-310 0 0 6.95259e-310 6.95259e-310 6.95259e-310 6.95259e-310 0 1.58856e-310 9.88131e-324 3.76832e-317 0 0 
      0 0 0 0 0 0 0 0 0 6.95259e-310 0 0 0 6.95259e-310 0 0 0 6.95259e-310 0 0 0 6.95259e-310 0 0 0 6.95259e-310 0 0 0 
      2.122e-314 0 0 0 6.95259e-310 0 0 6.95259e-310 6.95259e-310 6.95259e-310 6.95259e-310 0 1.58857e-310 1.99825e-316 
      3.76832e-317 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
      3.15399e-294 1.9986e-316 3.75609e-317 3.75602e-317 3.75602e-317 3.77014e-317 1.61895e-319
      ```