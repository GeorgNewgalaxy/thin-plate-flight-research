from lib import *
import matplotlib.pyplot as plt

reactants = ["G1", "G2", "G1*", "R", "R2", "P"]
reactions = [Reaction(R("G1"), R("G1*"), 0.5),
             Reaction(R("G1*"), R("G1"), 0.5),
             Reaction(R("G1*") + R("G2"), R("G1") + R("R"), 1),
             Reaction(R("G1") + R("R"), R("P") + R("R"), 1),
             Reaction(R("G1*") + R("R"), R("P") + R("R") + R("R"), 1),
             Reaction(R("G2") + R("R"), R("P") + R("R"), 0.1),
             Reaction(R("R"), R("R2"), 0.3)]

concentrations = dict(zip(reactants, [0]*len(reactants)))
concentrations["G1"] = 1

concentrations_with_time = dict(zip(reactants, [[]]*len(reactants)))

dt = 0.01
iterations = 8000
for i in range(0, iterations):
    if i == 1500:
        concentrations["G2"] = 0.1
    concentrations_new = concentrations.copy()
    for key in reactants:
        cwtk = concentrations_with_time[key].copy()
        cwtk.append(concentrations[key])
        concentrations_with_time[key] = cwtk.copy()
    for reagent in reactants:
        total_k = 0
        for reaction in reactions:
            c = reaction.reagents.reagents.count(reagent)
            if reaction.reagents.reagents.count(reagent) > 0:
                total_k -= reaction.calculate_rate(concentrations, reactants)*c
        for reaction in reactions:
            c = reaction.products.reagents.count(reagent)
            if reaction.products.reagents.count(reagent) > 0:
                total_k += reaction.calculate_rate(concentrations, reactants)*c
                two = 2

        concentrations_new[reagent] += dt*total_k
    concentrations = concentrations_new.copy()

for key in reactants:
    plt.plot(range(0, iterations), concentrations_with_time[key], label=key)

plt.legend()
plt.show()
