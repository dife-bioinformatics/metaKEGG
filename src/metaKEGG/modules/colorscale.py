from ..helpers.helpfunctions import get_colors_from_colorscale

colors_list = get_colors_from_colorscale(['tab20' , 'tab20b', 'tab20c',
                                          'Pastel1', 'Pastel2',
                                          'Paired', 'Accent', 'Dark2',
                                          'Set1', 'Set2', 'Set3'] , skip=2)
# colors_list.remove('#7f7f7f') # remove gray to avoid mix-up with background
colors_list.remove('#c7c7c7') # remove gray to avoid mix-up with background