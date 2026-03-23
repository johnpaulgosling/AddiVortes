linters <- lintr::linters_with_defaults(
  object_name_linter = NULL,
  object_length_linter = NULL,
  indentation_linter = NULL,
  commented_code_linter = NULL,
  return_linter = NULL,
  line_length_linter(120L)
)
