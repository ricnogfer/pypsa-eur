
# import necessary Python modules
import geopandas
import pandas
import logging
import yaml



def build_afforestation_potentials(config_yaml, network_geojson, nuts2_geojson, afforestation_corine_potentials_csv_file, afforestation_nuts2_rates_csv_file, output_csv_file, log):

    # load config yaml file representing the PyPSA-Eur configuration
    handle = open(config_yaml)
    config = yaml.safe_load(handle)
    handle.close()


    # configure log mechanism
    if log is True:
        logging.basicConfig(level = config["logging"]["level"])
        logger = logging.getLogger(__name__)
        logger.info("Calculate afforestation potentials")


    # load files
    network = geopandas.read_file(network_geojson)
    nuts2 = geopandas.read_file(nuts2_geojson)
    corine_potentials = pandas.read_csv(afforestation_corine_potentials_csv_file).set_index("node")
    nuts2_rates = pandas.read_csv(afforestation_nuts2_rates_csv_file).set_index("NUTS2")
    
    
    # harmonize NUTS2 indexes against the network indexes
    nuts2["NUTS_ID"] = nuts2["NUTS_ID"].apply(lambda x: "GR%s" % x[2:] if x[:2] == "EL" else x)
    nuts2["NUTS_ID"] = nuts2["NUTS_ID"].apply(lambda x: "GB%s" % x[2:] if x[:2] == "UK" else x)
    nuts2_rates = nuts2_rates.rename(index = lambda x: "GR%s" % x[2:] if x[:2] == "EL" else x)
    nuts2_rates = nuts2_rates.rename(index = lambda x: "GB%s" % x[2:] if x[:2] == "UK" else x)


    # create data frame to store afforestation potential for each node
    data_frame = pandas.DataFrame(columns = ["node", "potential [t/ha]"])


    # iterate through PyPSA-Eur network regions (nodes)
    for i in range(len(network)):

        # get node name and geometry
        region = network.iloc[i]
        node_name = region["name"]
        node_geometry = geopandas.GeoSeries(region["geometry"], crs = 3035)

        # iterate through NUTS2 codes
        node_afforestation_potential = 0
        for j in range(len(nuts2)):

            # get NUTS2 row        
            nuts2_row = nuts2.iloc[j]
            nuts2_name = nuts2_row["NUTS_ID"]

            # check that NUTS2 belong to the node's country
            if node_name[:2] != nuts2_name[:2]:
                continue

            # get NUTS2 geometry
            nuts2_geometry = geopandas.GeoSeries(nuts2_row["geometry"], crs = 3035)

            # calculate proportion of intersection between NUTS2 geometry and node geometry
            intersection = nuts2_geometry.intersection(node_geometry.iloc[0])
            proportion = float(intersection.area.iloc[0]) / float(node_geometry.area.iloc[0])

            # calculate afforestation potential based on intersection and aggregate this potential to the node's afforestation potential
            node_afforestation_potential += (corine_potentials.loc[node_name]["potential [sqkm]"] * 100) * nuts2_rates.loc[nuts2_name]["affo rate (t/ha/y)"] * proportion

        # add node afforestation potential into data frame
        if log is True:
            logger.info("Node '%s' has an afforestation potential of %d [t/ha]" % (node_name, node_afforestation_potential))
        data_frame.loc[len(data_frame)] = [node_name, node_afforestation_potential]


    # save afforestation potentials into CSV file
    if log is True:
        logger.info("Save afforestation potentials into CSV file '%s'" % output_csv_file)
    data_frame.set_index("node", inplace = True)
    data_frame.to_csv(output_csv_file)



if __name__ == "__main__":

    # build and save afforestation potentials into CSV file
    if "snakemake" in globals():
        build_afforestation_potentials("config/config.yaml", snakemake.params["network_geojson"], snakemake.params["nuts2_geojson"], snakemake.input["afforestation_corine_potentials_csv_file"], snakemake.input["afforestation_nuts2_rates_csv_file"], snakemake.output["csv_file"], True)
    else:
        build_afforestation_potentials("config.yaml", "regions_onshore_base_s_39.geojson", "NUTS_RG_03M_2013_4326_LEVL_2.geojson", "afforestation_corine_potentials_s_39.csv", "afforestation_nuts2.csv", "afforestation_potentials_s_39.csv", True)


