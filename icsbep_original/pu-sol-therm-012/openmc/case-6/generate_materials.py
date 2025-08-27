import openmc

mats = openmc.Materials()

mat = openmc.Material(1)
mat.name = "Plutonium nitrate solution (105.0 g/L)"
mat.set_density('sum')
mat.add_nuclide('Pu239', 1.96582e-04)
mat.add_nuclide('Pu240', 4.98109e-05)
mat.add_nuclide('Pu241', 1.47655e-05)
mat.add_nuclide('Pu242', 2.98268e-06)
mat.add_nuclide('Am241', 1.60010e-06)
mat.add_element('N', 2.51393e-03)
mat.add_element('O', 3.72502e-02)
mat.add_nuclide('H1', 6.07598e-02)
mat.add_element('Fe', 2.48014e-05)
mat.add_element('Cr', 7.97553e-06)
mat.add_element('Ni', 5.65269e-06)
mat.add_s_alpha_beta('c_H_in_H2O')
mats.append(mat)

mat = openmc.Material(2)
mat.name = "Air"
mat.set_density('sum')
mat.add_nuclide('O16', 1.0784e-05)
mat.add_nuclide('N14', 4.3090e-05)
mats.append(mat)

mat = openmc.Material(3)
mat.name = "Stainless steel"
mat.set_density('sum')
mat.add_element('Fe', 6.1344e-02)
mat.add_element('Cr', 1.6472e-02)
mat.add_element('Ni', 8.1050e-03)
mats.append(mat)

mat = openmc.Material(4)
mat.name = "Lucoflex"
mat.set_density('sum')
mat.add_element('C', 2.7365e-02)
mat.add_nuclide('H1', 4.1047e-02)
mat.add_element('Cl', 1.3682e-02)
mats.append(mat)

mat = openmc.Material(5)
mat.name = "Water"
mat.set_density('sum')
mat.add_nuclide('H1', 6.6688e-02)
mat.add_element('O', 3.3344e-02)
mat.add_s_alpha_beta('c_H_in_H2O')
mats.append(mat)

mats.export_to_xml()
