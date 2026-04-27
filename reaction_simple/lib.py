class R:
    def __init__(self, reagents):
        if type(reagents) == list:
            self.reagents = reagents
        else:
            self.reagents = [reagents]

    def __add__(self, other):
        reagents = self.reagents.copy()
        reagents.extend(other.reagents)
        return R(reagents)

    def __repr__(self):
        return str(self.reagents)


class Reaction:
    def __init__(self, reagents, products, k):
        self.reagents = reagents
        self.products = products
        self.k = k

    def calculate_rate(self, concentrations, keys):
        product = 1
        for key in keys:
            product *= (concentrations[key])**(self.reagents.reagents.count(key))
        return product*self.k

    def __repr__(self):
        out = ""
        for i in self.reagents.reagents:
            out += i
            out += "+"
        out += "->"
        for i in self.products.reagents:
            out += i
            out += "+"
        return out
