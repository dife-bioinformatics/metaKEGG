from ..helpers.helpfunctions import get_colors_from_colorscale

colors_list = get_colors_from_colorscale(['tab10' , 'tab20b', 'tab20c',
                                          'Pastel1', 'Pastel2',
                                          'Paired', 'Accent', 'Dark2',
                                          'Set1', 'Set2', 'Set3'])
colors_list.remove('#7f7f7f') # remove gray to avoid mix-up with background