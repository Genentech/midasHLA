# "4digit_C-allele_C1-2"
# Convert the file as follows:
# Reduce alleles to 4 digits resolution, remove duplicates
dict_name <- "4digit_C-allele_C1-2"
dict_path <- "MiDAS/inst/extdata/"
print(sprintf("checking %s dict", dict_name))
dict <-read.table(
  paste0(dict_path, "Match_", dict_name, ".txt"),
  header = TRUE,
  stringsAsFactors = FALSE
)
dict4 <- dict
dict4$allele <- reduceAlleleResolution(dict$allele, resolution = 4)
new_dict <- dict4[! duplicated(dict4$allele), ]
write.table(
  new_dict,
  file = paste0(dict_path, "Match_", dict_name, ".txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

# 4digit_B-allele_Bw
# Convert the file as follows:
# reduce all alleles to 4 digits, remove duplicates in case of conflicts favor Bw4 allele
# reduce all alleles to 6 digits, remove duplicates
# leave 8 digits alleles as is to preserve matchings for allele numbers with expression status indicated
# combine the 3 dictionaries to recreate full matching file
dict_name <- "4digit_C-allele_C1-2"
dict_path <- "MiDAS/inst/extdata/"
print(sprintf("checking %s dict", dict_name))
dict <-
  read.table(
    paste0(dict_path, "Match_", dict_name, ".txt"),
    header = TRUE,
    stringsAsFactors = FALSE
  )

dict4 <- dict
dict4$allele <- reduceAlleleResolution(dict4$allele, resolution = 4)
dict4 <- dict4[getAlleleResolution(dict4$allele) == 4, ] # Discard 8digit alleles
Bw_group <- vapply(
  X = unique(dict4$allele),
  FUN = function(x) {
    x <- unique(dict4[dict4$allele == x, "Bw"])
    ifelse(length(x) == 1, x, "Bw4")
  },
  FUN.VALUE = character(length = 1L)
)
dict4 <- data.frame(allele = unique(dict4$allele),
                    Bw = Bw_group,
                    stringsAsFactors = FALSE
)

dict6 <- dict
dict6$allele <- reduceAlleleResolution(dict6$allele, resolution = 6)
dict6 <- dict6[getAlleleResolution(dict6$allele) == 6, ] # Discard 8 & 4 digits alleles
dict6 <- dict6[! duplicated(dict6$allele), ]

dict8 <- dict
dict8 <- dict8[getAlleleResolution(dict8$allele) == 8, ]

new_dict <- bind_rows(dict4, dict6, dict8)
write.table(new_dict,
            file = paste0(dict_path, "Match_", dict_name, ".txt"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE
)

# 4digit_allele_Ggroup
dict_name <- "4digit_allele_Ggroup"
dict_path <- "MiDAS/inst/extdata/"
print(sprintf("checking %s dict", dict_name))
dict <-
  read.table(
    paste0(dict_path, "Match_", dict_name, ".txt"),
    header = TRUE,
    stringsAsFactors = FALSE
  )

print(sprintf("Percentage of duplicated entries: %.1f%%",
              100 * (1 - (length(unique(dict$allele)) / nrow(dict)))
))
duplicated_rows <- duplicated(dict$allele) | duplicated(dict$allele, fromLast = TRUE)
View(dict[duplicated_rows, ])

# alleles DRB1*01:01:01 and DRB1*01:16 are matched against more than one group,
# to solve this following matches are removed:
# DRB1*01:01:01 -> DRB5*01:01:01G
# DRB1*01:16 -> DRB5*01:01:01G
removed <- (dict$allele == "DRB1*01:01:01" | dict$allele == "DRB1*01:16") & dict$g_group == "DRB5*01:01:01G"
dict <- dict[! removed, ]

# DPA1 is wrongly named DPA, this is changed
dpa <- grepl("DPA", dict$allele)
dict$allele[dpa] <- gsub("^DPA", "DPA1", dict$allele[dpa])

write.table(dict,
            file = paste0(dict_path, "Match_", dict_name, ".txt"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE
)

#test
hla_calls <- readHlaCalls("test.txt")
hlaToVariable(hla_calls, dictionary = "4digit_allele_Ggroup")
