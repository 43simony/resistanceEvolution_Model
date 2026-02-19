library(httr)
library(jsonlite)
library(dplyr)
library(purrr)

citeList = read.csv2(file = header = T, sep = ',')





write.csv2(citeList, file = "allCites.csv")




get_citing_titles <- function(paper_id, limit = 1000) {
  
  # API endpoint: get citations for a paper
  base_url <- "https://api.semanticscholar.org/graph/v1/paper/"
  
  # We request: citing paper title + year + authors
  fields <- "title,year,authors"
  
  # Full URL
  url <- paste0(
    base_url, paper_id,
    "/citations?fields=", fields,
    "&limit=", limit
  )
  
  # Make request
  response <- GET(url)
  
  # Parse result
  data_raw <- fromJSON(content(response, "text", encoding = "UTF-8"))
  
  # Extract citing papers
  citing <- data_raw$data %>%
    map_df(~{
      tibble(
        title = .x$citingPaper$title,
        year  = .x$citingPaper$year
      )
    })
  
  return(citing)
}

# ---- Example usage ----
paper_id <- "8a1b5c4d2e973bba9e1c1faf8fbe44fabc123456"   # replace with your real ID
citing_titles <- get_citing_titles(paper_id)

print(citing_titles)