def remove_regex(str):
    regex_chars = "{}|()?¿*+|/.<>"
    for char in regex_chars:
        str = str.replace(char, "")
    return str
