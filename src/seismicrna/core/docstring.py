DEFAULT_DOCSTRING_INDENT = 2


def find_offset(docstring):
    docstring = clean_docstring(docstring)
    min_indent = 100
    for line in docstring.split('\n')[1:]:
        if is_line_empty(line):
            continue

        def count_indent(lin):
            for count, char in enumerate(lin):
                if char != '\t':
                    return count

        min_indent = min(min_indent, count_indent(line))
    return min_indent


def docstring_to_dict_google(docstring, offset=DEFAULT_DOCSTRING_INDENT):
    docstring = clean_docstring(docstring)

    out = {'header': docstring.split('\n')[0]}

    current_section = None
    for line in docstring.split('\n'):
        if is_line_empty(line):
            continue
        if line.startswith('\t' * (1 + offset)) and current_section is not None:
            out[current_section][line.strip().split(' ')[0]] = ' '.join(line.strip().split(' ')[1:])
        elif line.startswith('\t' * offset):
            current_section = line.strip()
            out[current_section] = {}

    return out


def dict_google_to_docstring(d, offset=DEFAULT_DOCSTRING_INDENT):
    docstring = ''
    if 'header' in d.keys():
        docstring = d.pop('header')

    docstring += '\n'
    for key in d.keys():
        docstring += '\n' + '\t' * offset + key + '\n'
        for item in d[key]:
            docstring += '\t' * (1 + offset) + item + ' ' + d[key][item] + '\n'

    return docstring


def style_child_takes_over_parent(prnt_doc, child_doc):
    docs = [prnt_doc, child_doc]
    for idx, doc in enumerate(docs):
        docs[idx] = docstring_to_dict_google(doc)
    prnt_doc, child_doc = docs

    for key in prnt_doc.keys():
        if key not in child_doc.keys():
            child_doc[key] = prnt_doc[key]
        else:
            if key == 'header':
                continue
            for item in prnt_doc[key]:
                if item not in child_doc[key].keys():
                    child_doc[key][item] = prnt_doc[key][item]

    return dict_google_to_docstring(child_doc)


def is_line_empty(line):
    return line.replace('\t', '').replace(' ', '') == ''


def clean_docstring(docstring):
    return docstring.replace('    ', '\t')
