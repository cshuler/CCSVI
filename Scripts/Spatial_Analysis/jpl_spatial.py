def func_zonal_stats(zone_fc, zone_field, raster_path, output_folder, stat_type):
    import os
    import arcpy
    import arcpy.sa as sa
    arcpy.CheckOutExtension("Spatial")
    arcpy.env.overwriteOutput = True

    # Sanitize raster base name for use in DBF field/table names (max 8 chars is safest for DBF)
    raster_base = os.path.splitext(os.path.basename(raster_path))[0]
    raster_base_clean = raster_base.replace(" ", "_")[:32]  # Max safe length for filename
    stat_type_clean = stat_type.lower()

    out_table_name = f"{raster_base_clean}_{stat_type_clean}.dbf"

    out_table_path = os.path.join(output_folder, "temp", out_table_name)
    
    out_temp_folder = os.path.join(output_folder, "temp")
    
    #Create OUTPUT folder which was pasased in as variable, and nested TEMP folder if they don't exist already
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(out_temp_folder, exist_ok=True)
    

    # Run zonal statistics
    sa.ZonalStatisticsAsTable(
        in_zone_data=zone_fc,
        zone_field=zone_field,
        in_value_raster=raster_path,
        out_table=out_table_path,
        statistics_type=stat_type,
        ignore_nodata="DATA"
    )
    
    return out_table_path  # return path to the table so you can print for confirmation


# #### Function to calculate whether the zone_fc vector layer has any contact with the vector_path layer

# In[4]:


def func_vector_match(zone_fc, zone_field, vector_path, output_folder, match_type):
    import os
    import arcpy
    arcpy.env.overwriteOutput = True

    # Sanitize base name
    vector_base = os.path.splitext(os.path.basename(vector_path))[0]
    vector_base_clean = vector_base.replace(" ", "_")[:32]
    match_type_clean = match_type.lower()

    # Output table name
    out_table_name = f"{vector_base_clean}_{match_type_clean}.dbf"
    out_table_path = os.path.join(output_folder, "temp", out_table_name)

    out_temp_folder = os.path.join(output_folder, "temp")
    
    #Create OUTPUT folder which was pasased in as variable, and nested TEMP folder if they don't exist already. Clean them if they do
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(out_temp_folder, exist_ok=True)

    # Create feature layers
    zone_layer = "zone_layer"
    vector_layer = "vector_layer"
    arcpy.MakeFeatureLayer_management(zone_fc, zone_layer)
    arcpy.MakeFeatureLayer_management(vector_path, vector_layer)

    # New code to dissolve
    dissolved_vector = os.path.join(out_temp_folder, "dissolved_temp.shp")
    # Delete existing output to avoid overwrite issues
    if arcpy.Exists(dissolved_vector):
        arcpy.management.Delete(dissolved_vector)
    arcpy.management.Dissolve(vector_path, dissolved_vector)
    
    # Perform spatial join
    temp_join = "in_memory/temp_join"
    arcpy.analysis.SpatialJoin(
        target_features=zone_layer,
        join_features=dissolved_vector,
        out_feature_class=temp_join,
        join_operation="JOIN_ONE_TO_ONE",
        join_type="KEEP_ALL",
        match_option=match_type.upper()  # INTERSECT, WITHIN, CONTAINS, etc.
    )

    # Add a binary 1/0 field indicating presence/absence of match
    fieldname_short = f"{match_type_clean}"[:10]  # Ensure DBF-safe
    arcpy.management.AddField(temp_join, fieldname_short, "SHORT")

    # Determine if match occurred by checking for nulls in join fields
    arcpy.management.CalculateField(
        in_table=temp_join,
        field=fieldname_short,
        expression="0 if !Join_Count! is None or !Join_Count! == 0 else 1",
        expression_type="PYTHON3"
    )

    # Export table
    arcpy.conversion.TableToTable(temp_join, out_temp_folder, out_table_name)

    return out_table_path


# #### Function to calculate % overlap between zone_fc vector layer has any contact with the vector_path layer

# In[5]:


def func_vector_pct(zone_fc, zone_field, vector_path, output_folder):
    import os
    import arcpy

    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")

    # Clean name for output (based on intersecting feature class)
    vector_base = os.path.splitext(os.path.basename(vector_path))[0]
    out_name = f"{vector_base[:32]}_pct.shp"
    out_path = os.path.join(output_folder, "temp", out_name)

    out_temp_folder = os.path.join(output_folder, "temp")
    
    #Create OUTPUT folder which was pasased in as variable, and nested TEMP folder if they don't exist already. Clean them if they do
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(out_temp_folder, exist_ok=True)

    # New code to dissolve
    dissolved_vector = os.path.join(out_temp_folder, "dissolved_temp.shp")
    # Delete existing output to avoid overwrite issues
    if arcpy.Exists(dissolved_vector):
        arcpy.management.Delete(dissolved_vector)
    arcpy.management.Dissolve(vector_path, dissolved_vector)
    
    # Step 1: Intersect zone with input vector (attributes from zone_fc retained)
    arcpy.analysis.Intersect([zone_fc, dissolved_vector], out_path, "ALL")

    # Step 2: Add overlap area (short field name)
    arcpy.management.AddField(out_path, "ov_area", "DOUBLE")
    arcpy.management.CalculateGeometryAttributes(out_path, [["ov_area", "AREA_GEODESIC"]])

    # Step 3: Get zone area from original features (via temporary copy)
    #zone_copy = os.path.join(output_folder, "zone_temp.shp")
    #arcpy.conversion.FeatureClassToFeatureClass(zone_fc, output_folder, out_name)
    zone_copy_name = "zone_temp.shp"
    zone_copy = os.path.join(out_temp_folder, zone_copy_name)
    arcpy.conversion.FeatureClassToFeatureClass(zone_fc, out_temp_folder, zone_copy_name)

    arcpy.management.AddField(zone_copy, "zn_area", "DOUBLE")
    arcpy.management.CalculateGeometryAttributes(zone_copy, [["zn_area", "AREA_GEODESIC"]])

    # Step 4: Join zone area to intersected features
    arcpy.management.JoinField(out_path, zone_field, zone_copy, zone_field, ["zn_area"])

    # Step 5: Add and calculate percent field
    arcpy.management.AddField(out_path, "pct", "FLOAT")
    arcpy.management.CalculateField(
        out_path, "pct",
        expression="(!ov_area! / !zn_area!) * 100 if !zn_area! else 0",
        expression_type="PYTHON3"
    )

    #Step 6: Delete the temporary zone shapefile
    arcpy.management.Delete(zone_copy)
    
    return out_path


# #### Function to calculate ACREAGE overlap between zone_fc vector layer has any contact with the vector_path layer

# In[6]:


def func_vector_acres(zone_fc, zone_field, vector_path, output_folder):
    import os
    import arcpy

    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")

    # Clean name for output shapefile
    vector_base = os.path.splitext(os.path.basename(vector_path))[0]
    out_name = f"{vector_base[:32]}_acres.shp"
    out_path = os.path.join(output_folder, "temp", out_name)

    out_temp_folder = os.path.join(output_folder, "temp")
    
    #Create OUTPUT folder which was pasased in as variable, and nested TEMP folder if they don't exist already
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(out_temp_folder, exist_ok=True)


    # Delete existing output to avoid overwrite issues
    if arcpy.Exists(out_path):
        arcpy.management.Delete(out_path)

    # New code to dissolve
    dissolved_vector = os.path.join(out_temp_folder, "dissolved_temp.shp")
    # Delete existing output to avoid overwrite issues
    if arcpy.Exists(dissolved_vector):
        arcpy.management.Delete(dissolved_vector)
    arcpy.management.Dissolve(vector_path, dissolved_vector)
    
    # Intersect zone features with vector_path (attributes from zone_fc retained)
    arcpy.analysis.Intersect([zone_fc, dissolved_vector], out_path, "ALL")

    # Add area field (in acres)
    arcpy.management.AddField(out_path, "acres", "DOUBLE")
    arcpy.management.CalculateGeometryAttributes(out_path, [["acres", "AREA_GEODESIC"]], area_unit="ACRES")

    return out_path


# In[7]:


def func_spatialcount(zone_fc, zone_field, vector_path, output_folder, radius):
    import os
    import arcpy

    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")

    # Clean name for output shapefile
    vector_base = os.path.splitext(os.path.basename(vector_path))[0]
    out_name = f"{vector_base[:32]}_count.shp"
    out_path = os.path.join(output_folder, "temp", out_name)

    out_temp_folder = os.path.join(output_folder, "temp")
    
    #Create OUTPUT folder which was pasased in as variable, and nested TEMP folder if they don't exist already. Clean them if they do
    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(out_temp_folder, exist_ok=True)

    # Delete existing output to avoid overwrite issues
    if arcpy.Exists(out_path):
        arcpy.management.Delete(out_path)

    # Step 1: Perform spatial join to count features within radius
    spatial_join_temp = os.path.join(output_folder, "temp_spjoin.shp")
    arcpy.analysis.SpatialJoin(
        target_features=zone_fc,
        join_features=vector_path,
        out_feature_class=spatial_join_temp,
        join_operation="JOIN_ONE_TO_MANY",
        join_type="KEEP_ALL",
        match_option="WITHIN_A_DISTANCE",
        search_radius=radius,
        distance_field_name=None
    )

    # Step 2: Dissolve to get count of joined features per zone
    arcpy.analysis.Statistics(
        in_table=spatial_join_temp,
        out_table=out_path,
        statistics_fields=[["JOIN_FID", "COUNT"]],
        case_field=zone_field
    )

    # Step 3 (optional): Delete temp join file
    arcpy.management.Delete(spatial_join_temp)

    return out_path


# #### Define a function to combine the multiple DBF files into one CSV. The field name in the CSV is taken from the DBF name, which originally was taken from the raster:

# In[8]:


def func_combine_tables(in_folder, out_folder, key_field, output_csv):
    import os
    import numpy as np
    import pandas as pd
    import geopandas as gpd
    import arcpy

    #key_field = "GEOIDFQ"
    #output_csv = "combined_zonal_stats.csv"
    #create list of DBF files in the in_folder
    dbf_files = [f for f in os.listdir(in_folder) if f.lower().endswith(".dbf")]
    merged_df = None

    for dbf in dbf_files:
        dbf_path = os.path.join(in_folder, dbf)
        
        # Get all non-geometry, non-OID fields
        fields = [f.name for f in arcpy.ListFields(dbf_path) if f.type not in ("Geometry", "OID")]

        if key_field not in fields or len(fields) < 2:
            continue  # Skip if missing key or only one useful field

        stat_field = fields[-1]  # Use rightmost field
        table = arcpy.da.TableToNumPyArray(dbf_path, [key_field, stat_field])
        df = pd.DataFrame(table)

        # Rename stat field to the DBF base name
        stat_name = os.path.splitext(dbf)[0]
        df = df.rename(columns={stat_field: stat_name})

        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on=key_field, how="outer")

    # Save combined CSV
    if merged_df is not None:
        output_csv_path = os.path.join(out_folder, output_csv)
        merged_df.to_csv(output_csv_path, index=False)
        return output_csv_path
    else:
        raise ValueError("No valid DBFs found or no data to merge.")

    return output_csv_path


# #### Define a function to join the CSV back onto the CBG shapefile as a new output:

# In[9]:


def func_csv_join_to_shp(shapefile_cbg_path, merged_csv):
    import os
    import arcpy

    arcpy.env.overwriteOutput = True

    # Output path based on CSV folder
    output_folder = os.path.dirname(merged_csv)
    output_shp = os.path.join(output_folder, "combined_zonal_shapes.shp")

    # Assume GEOIDFQ is the common key
    arcpy.management.MakeFeatureLayer(shapefile_cbg_path, "cbg_layer")
    arcpy.management.MakeTableView(merged_csv, "csv_table")

    arcpy.management.AddJoin("cbg_layer", "GEOIDFQ", "csv_table", "GEOIDFQ", "KEEP_COMMON")

    arcpy.management.CopyFeatures("cbg_layer", output_shp)

    # Clean up layers
    arcpy.management.Delete("cbg_layer")
    arcpy.management.Delete("csv_table")

    return output_shp


def func_clean_folder(output_folder_path):
    import os
    for filename in os.listdir(output_folder_path):
        file_path = os.path.join(output_folder_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                import shutil
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Could not delete {file_path}: {e}")

print("Johann's spatial functions loaded successfully!")            
