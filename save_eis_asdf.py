import sunpy.map
import matplotlib.pyplot as plt
import astropy.units as u
import asdf

eis_map =sunpy.map.Map("eis_20151018_102719_composition_sis.asdf")

eis_map.peek()


with asdf.open('eis)20151018_102719_composition_sis.asdf') as f:
  composition= f'['composition'] 
  print(composition)

