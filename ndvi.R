ndviTable <- read_csv('raw_data/ndviStats.csv', col_types = 'ccc') %>%
  mutate(across(everything(), .funs = str_sub(str_sub(., -1), end = -1)))
  separate(ndviMean, 'meanSpring', 'meanFall', sep = ',') %>%
  separate(ndviMedian, 'medianSpring', 'medianFall', sep = ' ') %>%
  separate(ndviStdDev, 'sdSpring', 'sdFall', sep = ' ') %>%
  mutate(across(.fns = ~ str_extract(., '[:digit:]\\.[:digit:]+')))


