class Rule:
    def __init__(self, reactants, products, k):
        self.reactants = reactants
        self.products = products
        self.k = k

    def __repr__(self):
        return f'|Reactants: {self.reactants}, Products: {self.products}, Coefficient: {self.k}|'


def read_rules(txt):
    lines = txt.split("\n")
    rules_array = []
    for i in range(0, len(lines)):
        reactants_line = []
        products_line = []
        k_line = []
        line = lines[i]

        string_append = False
        string_to_read = ""
        have_seen_curly_bracket = False

        for char_idx in range(0, len(line)):
            char = line[char_idx]
            if string_append:
                string_to_read += char
            if char == "[":
                string_to_read = ""
                string_append = True
            if char == "]":
                string_to_read = string_to_read[0:-1]
                if have_seen_curly_bracket:
                    products_line.append(string_to_read)
                else:
                    reactants_line.append(string_to_read)
                string_to_read = ""
                string_append = False

            if char == "{":
                have_seen_curly_bracket = True
                string_to_read = ""
                string_append = True
            if char == "}":
                string_to_read = string_to_read[0:-1]
                k_line.append(string_to_read)
                string_to_read = ""
                string_append = False
        rules_array.append(Rule(reactants_line, products_line, float(k_line[0])))
    return rules_array


def read_starting_conditions(txt):
    conditions = {}
    lines = txt.split("\n")
    for i in range(0, len(lines)):
        line = lines[i]
        name = ""
        expression = ""
        append_to_string = True
        for char_idx in range(0, len(line)):
            char = line[char_idx]
            if char == ":":
                append_to_string = False
            if append_to_string:
                name += char
            elif char != ":":
                expression += char
        conditions[name] = expression
    return conditions
