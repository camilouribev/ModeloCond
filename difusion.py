import pandas as pd
from itertools import combinations
import math


def difusion(comp, gases, P, T_g, Frac_n2, frac_comp, masa_total):
    comp_mezcla = comp.keys()

    LJ_mezcla = gases.molecular_diameters  # char length Lennard Jones[A]
    epsilon_mezcla = gases.Stockmayers  # potential depth over cte boltzmann[K]
    pesos_mol_mezcla = gases.MWs  # [g/mol]

    F_m = []  # fracción másica

    for compound in comp_mezcla:
        if compound == 'Nitrogen':
            F_m.append(Frac_n2)
        else:
            F_m.append(frac_comp)

    df = pd.DataFrame(list(
        zip(comp_mezcla, LJ_mezcla, epsilon_mezcla, pesos_mol_mezcla, F_m)),
        columns=['compound', 'σ [Å]', 'ϵ/k', '[g/mol]', 'Mass Fraction']
        ).dropna(axis=0)

    P_bar = P*1e-5

    A = 1.06036
    B = 0.15610
    C = 0.19300
    D = 0.47635
    E = 1.03587
    F = 1.52996
    G = 1.76474
    H = 3.89411

    diff_binary = []
    combi = list(zip(list(combinations(df['compound'], 2)),
            list(combinations(df['σ [Å]'], 2)),
            list(combinations(df['ϵ/k'], 2)),
            list(combinations(df['[g/mol]'], 2))))
    df.set_index('compound', inplace=True)
    df['Mass Flow [g]/h'] = df['Mass Fraction']*masa_total
    df['Moles /h'] = df['Mass Flow [g]/h']/df['[g/mol]']
    df['Molar Fraction'] = df['Moles /h']/df['Moles /h'].sum()

    for pair in combi:

        M_ab = (2/((1/pair[3][0])+(1/pair[3][1])))
        sigma = (pair[1][0]*pair[1][1])/2
        epsilon = (pair[2][0]*pair[2][1])**(1/2)

        T_ast = T_g/epsilon
        # Integral de colisión

        omega_d = (A/(T_ast**B))+(C/math.exp(D*T_ast))+(E/math.exp(F*T_ast))+(G/math.exp(H*T_ast))

        # Coeficiente de difusión [m²/s]

        Diff = 0.00266*(T_g**(3/2))/(P_bar*M_ab**(1/2)*sigma**2*omega_d)*0.0001

        diff_binary.append([pair[0][0], pair[0][1], Diff])

    diff_matrix = {}

    for chem in comp_mezcla:

        diff_matrix[chem] = [x+[df.loc[chem]['Mass Fraction']]for x in diff_binary if x[0] == chem or x[1] == chem]
        diff_matrix = {i: j for i, j in diff_matrix.items() if j != []}
        diff_in_mix = dict((key, []) for key in diff_matrix.keys())

    for k, v in diff_matrix.items():
        for pair in diff_matrix[k]:
            pair.remove(k)
            diff_in_mix[k].append(df.loc[pair[0]]['Molar Fraction']/pair[1])

    for k, v in diff_in_mix.items():
        diff_in_mix[k] = (1 - df.loc[k]['Mass Fraction'])/sum(diff_in_mix[k])
    df

    sub_mass = df[df['[g/mol]'] != df['[g/mol]']['Nitrogen']]
    M_v = sub_mass['[g/mol]'].mean()

    M_red_vap = (2/((1/M_v) + (1/df['[g/mol]']['Nitrogen'])))

    return (diff_in_mix, df, M_red_vap)
