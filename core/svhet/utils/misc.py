def is_between(pos, start, end):
    '''Check if pos is in [start, end]'''
    return pos >= start and pos <= end

def print_separator_line() -> None:
    import shutil
    terminal_width = shutil.get_terminal_size().columns
    print("-" * terminal_width)