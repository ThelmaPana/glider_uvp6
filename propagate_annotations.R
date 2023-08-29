#--------------------------------------------------------------------------#
# Project: glider_uvp6
# Script purpose: Propagate annotations (morphocluster + manual) from morphocluster project (6433) and science project (7587) to original projects.
# Date: 18/08/2023
# Author: Thelma Panaïotis
#--------------------------------------------------------------------------#

library(ecotaxarapi)
library(tidyverse)
library(lubridate)
library(reticulate)
pd <- import("pandas")


## Get project 6433 ----
#--------------------------------------------------------------------------#
# In this project, we want:
# - objects annotated in morphocluster
# - objects manually validated, overwritting morphocluster annotations

# Fields to extract
fields <- "obj.orig_id,fre.morphocluster_label,obj.classif_qual,obj.classif_when,obj.classif_who,txo.display_name,txo.id"

# Get objects
resp_6433 <- get_object_set(
  6433,
  ProjectFilters = ProjectFilters(statusfilter = "PV"),
  fields = fields,
  )

 # Rename details columns
colnames(resp_6433$details) <- str_split(fields, ",")[[1]]

# Convert to dataframe
df_6433 <- resp_6433 %>% 
  as.data.frame() %>% 
  as_tibble()

# Keep relevant columns and only validated objects
df_6433 <- df_6433 %>% 
  select(
    orig_id = details.obj.orig_id, 
    morphocluster_label = details.fre.morphocluster_label, 
    classif_qual = details.obj.classif_qual, 
    taxon = details.txo.display_name, 
    taxon_id = details.txo.id,
    classif_when = details.obj.classif_when,
    classif_who = details.obj.classif_who
    ) %>% 
  filter(classif_qual == "V") %>% 
  mutate(
    classif_when = ymd_hms(classif_when), # convert to datetime
    classif_qual = "validated" # reformat validation status
    ) 


# Flag morphocluster validations
df_6433 <- df_6433 %>% mutate(from_morphocluster = classif_when == ymd_hms("2022-08-17 17:49:45")) # Morphocluster validations at 2022-08-17 17:49:45


## Get project 7587 ----
#--------------------------------------------------------------------------#
# In this project, we want:
# - objects manually validated, overwritting morphocluster annotations (morphocluster annotations are identical to those in project 6433)

# Get objects
resp_7587 <- get_object_set(
  7587,
  ProjectFilters = ProjectFilters(statusfilter = "PV"),
  fields = fields,
)

# Rename details columns
colnames(resp_7587$details) <- str_split(fields, ",")[[1]]

# Convert to dataframe
df_7587 <- resp_7587 %>% 
  as.data.frame() %>% 
  as_tibble()

# Keep relevant columns and only validated objects
df_7587 <- df_7587 %>% 
  select(
    orig_id = details.obj.orig_id, 
    morphocluster_label = details.fre.morphocluster_label, 
    classif_qual = details.obj.classif_qual, 
    taxon = details.txo.display_name, 
    taxon_id = details.txo.id,
    classif_when = details.obj.classif_when,
    classif_who = details.obj.classif_who
  ) %>% 
  filter(classif_qual == "V") %>% 
  mutate(
    classif_when = ymd_hms(classif_when), # convert to datetime
    classif_qual = "validated" # reformat validation status
  ) 


# Flag morphocluster validations
df_7587 <- df_7587 %>% mutate(from_morphocluster = classif_when == ymd_hms("2022-08-17 17:49:45")) # Morphocluster validations at 2022-08-17 17:49:45


## Process morphocluster annotations ----
#--------------------------------------------------------------------------#
# Flag morphocluster annotations and clean label
df_morpho <- df_6433 %>% 
  filter(from_morphocluster) %>% 
  mutate(
    classif_who = "",
    classif_mail = "",
    from_morphocluster = as.character(from_morphocluster),
    morphocluster_label = str_remove(morphocluster_label, "/80")
    )


## Process manual validations ----
#--------------------------------------------------------------------------#
# Make sure not to have duplicated annotations from 7584 and 6433
# Set the morphocluster annotation and flag to ""
df_man <- df_6433 %>% 
  filter(!from_morphocluster) %>% 
  filter(!orig_id %in% df_7587$orig_id) %>% 
  bind_rows(df_7587 %>% filter(!from_morphocluster)) %>% 
  mutate(
    classif_who = 469,
    morphocluster_label = "",
    from_morphocluster = "",
    classif_who = "Thelma Panaiotis",
    classif_mail = "thelma.panaiotis@imev-mer.fr"
  )


## Get internal object ids from original projects ----
#--------------------------------------------------------------------------#
# We need internal objects IDs from original projects to update objects
# Query relevant projects with the API

# Fields to extract
fields <- "obj.orig_id,obj.classif_qual,obj.classif_when,obj.classif_who,txo.display_name"

# Projects to extract
projects <- c(4028, 4041, 4065, 4171, 4211, 4282, 4360, 4448, 4545, 4551)

# Initiate empty tibble
df_origs <- tibble()

# Loop over projects to extract them
for (proj in projects){
  resp <- get_object_set(
    proj,
    ProjectFilters = ProjectFilters(statusfilter = "PV"),
    fields = fields,
  )
  
  # Rename details columns
  colnames(resp$details) <- str_split(fields, ",")[[1]]
  
  # Convert to dataframe
  df <- resp %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    select(object_id = object_ids, orig_id = details.obj.orig_id, project_ids)
  
  # Store with other objects
  df_origs <- bind_rows(df_origs, df)
  
}


## Join internal object IDs from original projects with new annotations ----
#--------------------------------------------------------------------------#
df <- df_morpho %>% 
  bind_rows(df_man) %>% 
  left_join(df_origs) %>% # join with internal ids
  select( # keep relevant columns
    orig_id,
    project_id = project_ids,
    object_id,
    object_annotation_category = taxon,
    object_annotation_person_name = classif_who,
    object_annotation_person_email = classif_mail,
    object_annotation_status = classif_qual,
    object_from_morphocluster = from_morphocluster,
    object_morphocluster_label = morphocluster_label
  ) 


## Ecotaxa formatting ----
#--------------------------------------------------------------------------#
# Prepare dataframe with ecotaxa formatting

#df <- df %>% select(-c(orig_id, project_id))
#
## first_row as Dataframe row with appropriate headers
#first_row = c('[t]', '[t]', '[t]', '[t]', '[t]', '[t]', '[t]')
#first_row = t(pd$DataFrame(first_row))
#colnames(first_row) <- colnames(df)
#
## concat first_row and dataframe
#df <- rbind(first_row, df)
#
## save table
#write_tsv(df, 'data/annotation_propagation/ecotaxa_manual_and_morphocluster.tsv')

# This does not work, it creates the new columns but they remain empty
# Let’s do the updates from R with the API


## Get list of taxa and internal IDs to update through the API ----
#--------------------------------------------------------------------------#
taxa <- df_6433 %>% 
  bind_rows(df_7587) %>% 
  select(object_annotation_category = taxon, object_annotation_category_id = taxon_id) %>% 
  unique()


## Update classifications and metadata using API ----
#--------------------------------------------------------------------------#
# Loop over projects
for (proj in projects){
  print(proj)
    
  # Get objects for this project
  up_proj <- df %>% filter(project_id == proj) %>% select(-c(orig_id, project_id))

  # Get Morphocluster annotations
  up_morpho <- up_proj %>% filter(object_from_morphocluster == "TRUE")
  print("morpho")
  
  # Loop over annotations and update objects
  for (taxon in unique(up_morpho$object_annotation_category)){
    up_morpho_tax <- up_morpho %>% 
      filter(object_annotation_category == taxon) %>% 
      left_join(taxa, by = join_by(object_annotation_category))

    # Update taxonomy 
    classify_object_set(
      ClassifyReq(
        target_ids = up_morpho_tax %>% pull(object_id),
        classifications = up_morpho_tax %>% pull(object_annotation_category_id) %>% as.numeric(),
        wanted_qualification = "V"
      ))
    
    # Update metadata
    update_object_set(
      BulkUpdateReq(
        target_ids = up_morpho_tax %>% pull(object_id),
        updates = list(
          ColUpdate(ucol = "from_morphocluster", uval = unique(up_morpho_tax$object_from_morphocluster)), 
          ColUpdate(ucol = "morphocluster_label", uval = unique(up_morpho_tax$object_morphocluster_label))
        )))
  }
  
  # Get manual annotations
  up_man <- up_proj %>% filter(object_from_morphocluster == "")
  print("manual")
  
  # Loop over annotations and update taxonomy
  # No need to update metadata, annotation person is automatically updated
  for (taxon in unique(up_man$object_annotation_category)){
    up_man_tax <- up_man %>% 
      filter(object_annotation_category == taxon) %>% 
      left_join(taxa, by = join_by(object_annotation_category))
    
    # Prepare taxomomy update
    if (nrow(up_man_tax) == 1){ # If there is only one object to update, pass an explicit list
      classify_req <- ClassifyReq(
        target_ids = list(up_man_tax %>% pull(object_id)),
        classifications = list(up_man_tax %>% pull(object_annotation_category_id) %>% as.numeric()),
        wanted_qualification = "V"
      )
    } else {
      classify_req <- ClassifyReq(
        target_ids = c(up_man_tax %>% pull(object_id)),
        classifications = up_man_tax %>% pull(object_annotation_category_id) %>% as.numeric(),
        wanted_qualification = "V"
      )
    }
    
    # Update taxonomy 
    classify_object_set(classify_req)
  }
}

