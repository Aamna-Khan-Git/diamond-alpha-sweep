"""
port_diamond.py - convert the authors' DIAMOnD.py (Python 2) to Python 3.

The original repository has no licence file, so instead of copying the code
into this repository we download it and apply a small set of mechanical
edits. None of the edits change the algorithm; they only change syntax:

  1. `import cPickle`            -> `import pickle as cPickle`  (renamed module)
  2. `dict.iteritems()`          -> `dict.items()`              (renamed method)
  3. `scipy.infty`               -> `np.inf`                    (removed from SciPy)
  4. `print x` / `print>>f, x`   -> `print(x)` / `print(x, file=f)`
  5. `lgn1 - [lgnk1 + lgk1]`     -> `lgn1 - (lgnk1 + lgk1)`
     The original subtracts a one-element *list*, which turns every p-value
     into a one-element NumPy array. The number is identical; using
     parentheses just keeps it a plain float.
  6. A backslash in one docstring ("\\sum") is escaped, to avoid a warning.
  7. `import scipy.special` is added because the code calls
     scipy.special.gammaln, and newer SciPy versions do not import that
     submodule automatically.

tests/test_port.py checks that the ported p-value equals SciPy's own
hypergeometric tail probability.
"""
import re

# matches a single- or double-quoted string literal
_STRING = re.compile(r"'(?:\\.|[^'\\])*'|\"(?:\\.|[^\"\\])*\"")


def _bracket_depth(text):
    """Count '(' minus ')' while ignoring brackets inside string literals."""
    code = _STRING.sub("", text)
    return code.count("(") - code.count(")")


def _convert_print_statements(lines):
    """Turn Python 2 print statements into print() calls.

    A print statement can continue over several lines when it has an open
    bracket, so we keep collecting lines until the brackets balance and then
    add the closing ')' at the end of the last line.
    """
    out, i = [], 0
    while i < len(lines):
        line = lines[i]
        match = re.match(r"^(\s*)print\s*(>>\s*(\w+)\s*,)?\s*(.*)$", line)
        is_statement = (
            match
            and not line.lstrip().startswith("print(")
            and not line.lstrip().startswith("print_")
        )
        if not is_statement:
            out.append(line)
            i += 1
            continue

        indent, _, target, body = match.groups()
        # collect continuation lines until brackets are balanced
        block = [body]
        depth = _bracket_depth(body)
        while depth > 0:
            i += 1
            block.append(lines[i].rstrip("\n"))
            depth += _bracket_depth(lines[i])
        text = "\n".join(b.rstrip("\n") for b in block)
        file_arg = f", file={target}" if target else ""
        out.append(f"{indent}print({text}{file_arg})\n")
        i += 1
    return out


def port_to_python3(source: str) -> str:
    """Return the Python 3 version of the original DIAMOnD.py source."""
    source = source.replace("import cPickle", "import pickle as cPickle")
    source = source.replace(
        "import scipy.stats", "import scipy.stats\nimport scipy.special"
    )
    source = source.replace(".iteritems()", ".items()")
    source = source.replace("scipy.infty", "np.inf")
    source = source.replace("lgn1 - [lgnk1 + lgk1]", "lgn1 - (lgnk1 + lgk1)")
    # a docstring contains "\sum", an invalid escape in Python 3.12+ (warning only)
    source = source.replace("p-val = \\sum", "p-val = \\\\sum")
    lines = source.splitlines(keepends=True)
    return "".join(_convert_print_statements(lines))
