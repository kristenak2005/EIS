import configparser

config = configparser.ConfigParser()

config.add_section('directories')
config.set('directories', 'data_dir', '/mnt/scratch/data/spruksk2/data/')
config.set('directories', 'main_dir', '/mnt/scratch/data/spruksk2/python_output')

with open(r"configfile.ini", 'w') as configfile:
    config.write(configfile)
