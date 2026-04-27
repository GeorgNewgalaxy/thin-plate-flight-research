import numpy as np
import matplotlib.pyplot as plt
from rule_processing import *

rules_raw_unread = open("rules", "r")
rules_raw = rules_raw_unread.read()
rules_raw_unread.close()


starting_conditions_raw_unread = open("starting_conditions", "r")
starting_conditions_raw = starting_conditions_raw_unread.read()
starting_conditions_raw_unread.close()

iterations = 10000

n_x = 100
n_t = 10

dx = 1/n_x
dt = 1/n_t

D = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001]

rules = read_rules(rules_raw)
starting_conditions = read_starting_conditions(starting_conditions_raw)
all_chemicals = list(starting_conditions.keys())

Vals = [[eval(starting_conditions[list(starting_conditions.keys())[i]])] for i in range(0, len(starting_conditions.keys()))]

for i in range(0, len(Vals)):
    assert D[i] * dt / dx**2 <= 0.5, f"Reactant {i} convergence failed"

for t in range(0, iterations):
    lasts = [[Vals[j][-1][i] for i in range(0, len(Vals[j][-1]))] for j in range(0, len(Vals))]

    news = [[0 for i in range(0, len(lasts[j]))] for j in range(0, len(lasts))]

    for i in range(1, n_x - 1):
        for j in range(0, len(Vals)):
            rules_with_j_reactant = []
            rules_with_j_product = []
            for rule in rules:
                for REPEAT in range(0, rule.reactants.count(all_chemicals[j])):
                    rules_with_j_reactant.append(rule)
                for REPEAT in range(0, rule.products.count(all_chemicals[j])):
                    rules_with_j_product.append(rule)

            product_of_concentration_reactant_list = []
            product_of_concentration_product_list = []
            for rule in rules_with_j_reactant:
                prod = 1
                for reactant in rule.reactants:
                    prod *= lasts[all_chemicals.index(reactant)][i]
                product_of_concentration_reactant_list.append(prod)

            for rule in rules_with_j_product:
                prod = 1
                for reactant in rule.reactants:
                    prod *= lasts[all_chemicals.index(reactant)][i]
                product_of_concentration_product_list.append(prod)

            s = 0
            for concentration_idx in range(0, len(rules_with_j_reactant)):
                s -= rules_with_j_reactant[concentration_idx].k*product_of_concentration_reactant_list[concentration_idx]
            for concentration_idx in range(0, len(rules_with_j_product)):
                s += rules_with_j_product[concentration_idx].k*product_of_concentration_product_list[concentration_idx]

            news[j][i] = lasts[j][i] + dt*D[j]*(lasts[j][i + 1] - 2*lasts[j][i] + lasts[j][i - 1])/(dx**2) + dt*s

    for i in range(0, len(Vals)):
        news[i][0] = news[i][1]
        news[i][n_x - 1] = news[i][n_x - 2]

    for i in range(0, len(Vals)):
        Vals[i].append(news[i])


ax = plt.axes(projection="3d")
x_data = np.arange(0, iterations, 1)
y_data = np.arange(0, n_x, 1)
X, Y = np.meshgrid(x_data, y_data)
Zs = [np.array([[Vals[i][t][x] for t in range(0, len(x_data))] for x in range(0, len(y_data))]) for i in range(0, len(Vals))]

chemical = "Y"
ax.plot_surface(X, Y, Zs[all_chemicals.index(chemical)])
plt.show()
