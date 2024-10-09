
# import necessary Python modules
import matplotlib.pyplot as plt
import geopandas
import pandas
import atlite
import logging
from atlite.gis import ExclusionContainer
from atlite.gis import shape_availability
from rasterio.plot import show



if __name__ == "__main__":

    # mock Snakemake in case Python module is called from a terminal
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_potentials")



    # get component to calculate the potential per node
    component = snakemake.params["component"]



    # configure log mechanism
    logging.basicConfig(level = snakemake.config["logging"]["level"])
    logger = logging.getLogger(__name__)
    logger.info("Calculate %s potentials" % component)



    # load geojson file representing the PyPSA-Eur network
    nodes_geojson = geopandas.read_file(snakemake.input["network_geojson"]).set_index("name")



    # select codes (specified in the "config.yaml" file) from Corine Land Cover (CLC) dataset
    excluder = ExclusionContainer(crs = 3035)
    excluder.add_raster(snakemake.input["corine_dataset"], codes = snakemake.config[component]["corine"], invert = True, crs = 3035)
    cell_area = excluder.res**2



    # calculate potential per node
    df = pandas.DataFrame(columns = ["node", "potential (t)"])
    for node in nodes_geojson.index:
        shape = nodes_geojson.to_crs(excluder.crs).loc[[node]].geometry
        band, transform = shape_availability(shape, excluder)
        selected_cells = band.sum()
        potential = selected_cells * cell_area / 1e6 * snakemake.config[component]["potential_per_sqkm"]   # in tonnes
        df.loc[len(df)] = [node, potential]
        #print("Node=%s" % node)
        #print("Area (km2)=%.2f" % (shape.geometry.area.sum() / 1e6))
        #print("Potential (t)=%.2f" % potential)
        #print()



    # save potentials into a CSV file
    logger.info("Save %s potentials into CSV file" % component)
    df.set_index("node", inplace = True)
    df.to_csv(snakemake.output["csv_file"])



    # save potentials into a PNG file
    logger.info("Save %s potentials into PNG file" % component)
    shape = nodes_geojson.to_crs(excluder.crs).geometry
    band, transform = shape_availability(shape, excluder)
    fig, ax = plt.subplots(figsize = (20, 23))
    ax.set_axis_off()
    shape.plot(ax = ax, color = "none")
    show(band, transform = transform, cmap = "Greens", ax = ax)
    plt.savefig(snakemake.output["png_file"])


