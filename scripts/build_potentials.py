
# import necessary Python modules
from atlite.gis import ExclusionContainer
from atlite.gis import shape_availability
from rasterio.plot import show
import matplotlib.pyplot as plt
import geopandas
import pandas
import atlite
import logging
import yaml



def build_potentials(config_yaml, network_geojson, corine_dataset, component, csv_file, png_file, log):

    # load config yaml file representing the PyPSA-Eur configuration
    handle = open(config_yaml)
    config = yaml.safe_load(handle)
    handle.close()


    # configure log mechanism
    if log is True:
        logging.basicConfig(level = config["logging"]["level"])
        logger = logging.getLogger(__name__)
        logger.info("Calculate %s potentials" % component)


    # load network geojson file representing the PyPSA-Eur network
    nodes_geojson = geopandas.read_file(network_geojson).set_index("name")


    # select CORINE Land Cover (CLC) codes (specified in the config yaml file)
    excluder = ExclusionContainer(crs = 3035)
    excluder.add_raster(corine_dataset, codes = config[component]["corine"], invert = True, crs = 3035)
    cell_area = excluder.res**2


    # calculate potential per node
    df = pandas.DataFrame(columns = ["node", "potential"])
    for node in nodes_geojson.index:
        shape = nodes_geojson.to_crs(excluder.crs).loc[[node]].geometry
        band, transform = shape_availability(shape, excluder)
        selected_cells = band.sum() * cell_area / 1e6   # in sqkm
        potential = selected_cells * config[component]["potential_per_sqkm"]
        df.loc[len(df)] = [node, potential]
        #print("Node=%s" % node)
        #print("Area (km2)=%.2f" % (shape.geometry.area.sum() / 1e6))
        #print("Potential=%.2f" % potential)
        #print()


    # save potentials into a CSV file
    if csv_file is not None:
        if log is True:
            logger.info("Save %s potentials into CSV file" % component)
        df.set_index("node", inplace = True)
        df.to_csv(csv_file)


    # save potentials into a PNG file
    if png_file is not None:
        if log is True:
            logger.info("Save %s potentials into PNG file" % component)
        shape = nodes_geojson.to_crs(excluder.crs).geometry
        band, transform = shape_availability(shape, excluder)
        fig, ax = plt.subplots(figsize = (20, 23))
        ax.set_axis_off()
        shape.plot(ax = ax, color = "none")
        show(band, transform = transform, cmap = "Greens", ax = ax)
        plt.savefig(png_file)



if __name__ == "__main__":


    # mock Snakemake in case Python module is called from a terminal
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_potentials")


    # build and save potentials into CSV and PNG files
    build_potentials("config.yaml", snakemake.input["network_geojson"], snakemake.input["corine_dataset"], snakemake.params["component"], snakemake.output["csv_file"], snakemake.output["png_file"], True)

