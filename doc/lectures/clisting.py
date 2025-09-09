clisting_filename = ""
clisting_cur_lineno = -1
clisting_cur_line = ""
clisting_all_lines = []


def clisting_read_next():
    global clisting_cur_line
    global clisting_cur_lineno

    clisting_cur_lineno += 1
    clisting_cur_line = clisting_all_lines[clisting_cur_lineno]


def clisting_read_prev():
    global clisting_cur_line
    global clisting_cur_lineno

    clisting_cur_lineno -= 1
    clisting_cur_line = clisting_all_lines[clisting_cur_lineno]


def clisting_read(line):
    while(True):
        try:
            clisting_read_next()
        except IndexError:
            raise Exception(f"Failed to find \"{line}\" in file {clisting_filename}")
        if clisting_cur_line.find(line) >= 0:
            break
    else:
        raise Exception(f"Line not found {line}")


def minted_output(start, end):
    while clisting_all_lines[start].strip() == "":
        start += 1
    while clisting_all_lines[end].strip() == "":
        end -= 1

    s1 = f"\\inputminted[firstline={start+1}, lastline={end+1}]"
    out = s1 + "{c++}{" + clisting_filename + "}"
    print(out)


def current_indent():
    l1 = len(clisting_cur_line)
    l2 = len(clisting_cur_line.lstrip())
    return clisting_cur_line[:l1-l2]


def clisting_open(fn):
    global clisting_all_lines
    global clisting_filename
    global clisting_cur_lineno
    global clisting_cur_line

    clisting_filename = fn
    fd = open(clisting_filename, 'r')
    clisting_all_lines = fd.readlines()
    clisting_cur_lineno = -1
    clisting_cur_line = ""
    fd.close()


def call_clisting_open(fn):
    clisting_open("../../src/" + fn)


def call_clisting_line(line):
    clisting_read(line)
    minted_output(clisting_cur_lineno, clisting_cur_lineno)


def call_clisting_block(line):
    clisting_read(line)
    l1 = clisting_cur_lineno
    indent = current_indent()
    needed_line = indent + "}"
    clisting_read_next()
    while not clisting_cur_line.startswith(needed_line):
        clisting_read_next()
    # fd = open("tmp.tmp", 'w')
    # fd.write(needed_line)
    # fd.write('\n-----\n')
    # fd.write(clisting_cur_line)
    # fd.close()
    l2 = clisting_cur_lineno
    minted_output(l1, l2)


def call_clisting_lines_range(line1, line2):
    clisting_read(line1)
    l1 = clisting_cur_lineno
    clisting_read(line2)
    l2 = clisting_cur_lineno
    minted_output(l1, l2)


def call_clisting_pass(line):
    clisting_read(line)


def call_clisting_to_start():
    clisting_open(clisting_filename)


def call_clisting_until(line):
    clisting_read_next()
    l1 = clisting_cur_lineno
    if line in clisting_cur_line:
        l2 = l1
    else:
        clisting_read(line)
        l2 = clisting_cur_lineno
    minted_output(l1, l2)


def call_clisting_before(line):
    clisting_read_next()
    l1 = clisting_cur_lineno
    clisting_read(line)
    clisting_read_prev()
    l2 = clisting_cur_lineno
    minted_output(l1, l2)


def clisting_func(what, arg1=0, arg2=0):
    if what == "open":
        call_clisting_open(arg1)
    elif what == "line":
        call_clisting_line(arg1)
    elif what == "lines-range":
        call_clisting_lines_range(arg1, arg2)
    elif what == "pass":
        call_clisting_pass(arg1)
    elif what == "to-start":
        call_clisting_to_start()
    elif what == "until":
        call_clisting_until(arg1)
    elif what == "until-close":
        call_clisting_until("}")
    elif what == "block":
        call_clisting_block(arg1)
    elif what == "before":
        call_clisting_before(arg1)
    else:
        raise Exception(f"Unknown clisting functional {what}")
