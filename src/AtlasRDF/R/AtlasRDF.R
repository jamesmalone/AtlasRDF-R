#library(SPARQL)
#library(hash)


###################
#functions to perform querying of Atlas data
###################


#######
#get all ensembl genes for an efo term for any species
#######
getAllEnsemblGenesForExFactor <- function(exfactor, limit = 0, endpoint="http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    
    limitC = ""
    if (limit != 0)
    if (is.numeric(limit) && ! is.null( grep ("/.",limit)))
    limitC = paste( " limit " , limit)             
    else
    warning ("limit should be an integer, limit omitted from the query")
    
    
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX efo:<http://www.ebi.ac.uk/efo/> \n",
            "PREFIX obo:<http://purl.obolibrary.org/obo/> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX sio:<http://semanticscience.org/resource/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            
            "SELECT distinct ?dbXref ?geneName ?ensemblid ?propertyValue \n",      
            "WHERE { \n",           
                "?expUri atlasterms:hasAnalysis ?analysis . \n",    
                "?analysis atlasterms:hasExpressionValue ?value . \n",       
                "?value atlasterms:hasFactorValue ?factor . \n",      
                "?value atlasterms:isMeasurementOf ?probe . \n",     
                "?probe atlasterms:dbXref ?dbXref . \n",
                "?dbXref rdf:type <http://rdf.ebi.ac.uk/terms/atlas/EnsemblDatabaseReference> . \n",
                "?dbXref rdfs:label ?geneName . \n",  
                "?dbXref dcterms:identifier ?ensemblid . \n",    
                "?factor atlasterms:propertyType ?propertyType . \n",       
                "?factor atlasterms:propertyValue ?propertyValue . \n",
                "?factor rdf:type " , exfactor , " .  \n",     
                "}  \n",
            limitC )
    
    message("Performing query please wait...")
    
    res <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("an error occured when trying to query for ensembl genes ", err)
            })#end tryCatch
    
    #res<-SPARQL(url=endpoint,query)
    return (res$results)    
}


##########
#get all ensembl genes for an efo term for a specified species only
##########
getSpeciesSpecificEnsemblGenesForExFactor <- function(exfactor, taxon, limit = 0, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    limitC = ""
    if (limit != 0)
    if (is.numeric(limit) && ! is.null( grep ("/.",limit)))
    limitC = paste( " limit " , limit)             
    else
    warning ("limit should be an integer, limit omitted from the query")
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX efo:<http://www.ebi.ac.uk/efo/> \n",
            "PREFIX obo:<http://purl.obolibrary.org/obo/> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX sio:<http://semanticscience.org/resource/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            
            "SELECT DISTINCT ?dbXref ?genename ?ensemblid \n",      
            "WHERE { \n",           
                "?expUri atlasterms:hasAnalysis ?analysis . \n",    
                "?analysis atlasterms:hasExpressionValue ?value . \n",       
                "?value atlasterms:hasFactorValue ?factor . \n",      
                "?value atlasterms:isMeasurementOf ?probe . \n",     
                "?probe atlasterms:dbXref ?dbXref . \n",
                "?dbXref rdf:type <http://rdf.ebi.ac.uk/terms/atlas/EnsemblDatabaseReference> . \n",
                "?dbXref rdfs:label ?genename . \n",  
                "?dbXref dcterms:identifier ?ensemblid . \n",
                "?dbXref atlasterms:taxon" , taxon , ". \n",    
                "?factor atlasterms:propertyType ?propertyType . \n",       
                "?factor atlasterms:propertyValue ?propertyValue . \n",
                "?factor rdf:type" , exfactor , " .  \n",     
                "}  \n",
            limitC )
    
    message("Performing query please wait...")
    
    res <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("an error occured when trying to query for species specific ensembl genes ", err)
            })#end tryCatch

    return (res$results)    
}



########
#get experiments where sample description contain specified term
########
getExperimentsByDescription <- function(searchterm, limit = 0, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    limitC = ""
    if (limit != 0)
    if (is.numeric(limit) && ! is.null( grep ("/.",limit)))
    limitC = paste( " limit " , limit)             
    else
    warning ("limit should be an integer, limit omitted from the query")
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX efo:<http://www.ebi.ac.uk/efo/> \n",
            "PREFIX obo:<http://purl.obolibrary.org/obo/> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX sio:<http://semanticscience.org/resource/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            
            "SELECT DISTINCT ?experiment WHERE \n", 
            "{?experiment a atlasterms:Experiment ; dcterms:description ?description ; \n", 
                "atlasterms:hasAssay \n", 
                "[atlasterms:hasSample \n", 
                        "[atlasterms:hasSampleCharacteristic \n", 
                                "[ atlasterms:propertyType ?propertyType ; atlasterms:propertyValue ?propertyValue] ] ] \n", 
                "filter regex (?description, \"",searchterm,"\") \n", 
                "}\n", 
            limitC, sep="")
    
    message("Performing query please wait...")

    res <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("an error occured when trying to query for experiment by description ", err)
            })#end tryCatch

    return (res$results)    
}


########
#get all genes for given experiment, expeirment should be sent in form of experiment ID e.g. E-GEOD-1085
########
getGenesForExperimentID <- function(experiment, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    experimenturi <- paste("atlas:",experiment, sep="")
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX efo:<http://www.ebi.ac.uk/efo/> \n",
            "PREFIX obo:<http://purl.obolibrary.org/obo/> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX sio:<http://semanticscience.org/resource/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            
            "SELECT distinct ?expressionValue ?dbXref ?pvalue ?propertyValue \n",       
            "WHERE {\n",           
                experimenturi,"atlasterms:hasAnalysis ?analysis . \n",    
                "?analysis atlasterms:hasExpressionValue ?value . \n",           
                "?value rdfs:label ?expressionValue . \n",     
                "?value atlasterms:pValue ?pvalue . \n",      
                "?value atlasterms:hasFactorValue ?factor . \n",      
                "?value atlasterms:isMeasurementOf ?probe . \n",     
                "?probe atlasterms:dbXref ?dbXref . \n",       
                "?factor atlasterms:propertyType ?propertyType . \n",       
                "?factor atlasterms:propertyValue ?propertyValue . \n",        
                "}")
    
    message("Performing query please wait...")

    res <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("an error occured when trying to query for genes for an experiment ID ", err)
            })#end tryCatch
    
    return (res$results)          
}



########
#get all genes for given experiment, expeirment should be sent in form of experiment URI e.g. <http://rdf.ebi.ac.uk/resource/atlas/E-GEOD-13396>
########
getGenesForExperimentURI <- function(experiment, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX efo:<http://www.ebi.ac.uk/efo/> \n",
            "PREFIX obo:<http://purl.obolibrary.org/obo/> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX sio:<http://semanticscience.org/resource/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            
            "SELECT distinct ?expressionValue ?dbXref ?pvalue ?propertyValue \n",       
            "WHERE {\n",           
            experiment,"atlasterms:hasAnalysis ?analysis . \n",    
            "?analysis atlasterms:hasExpressionValue ?value . \n",           
            "?value rdfs:label ?expressionValue . \n",     
            "?value atlasterms:pValue ?pvalue . \n",      
            "?value atlasterms:hasFactorValue ?factor . \n",      
            "?value atlasterms:isMeasurementOf ?probe . \n",     
            "?probe atlasterms:dbXref ?dbXref . \n",       
            "?factor atlasterms:propertyType ?propertyType . \n",       
            "?factor atlasterms:propertyValue ?propertyValue . \n",        
            "}")
    
    message("Performing query please wait...")
    
    res <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("an error occured when trying to query for genes for an experiment URI ", err)
            })#end tryCatch
    
    return (res$results)          
}


######
#query for experiments that contain a specific ensembl gene ID and which is reported as diff expressed
######
getExperimentURIsForGeneId <- function(geneid, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){    
    ensembluri <- paste("identifiers:", geneid, sep="")
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX efo:<http://www.ebi.ac.uk/efo/> \n",
            "PREFIX obo:<http://purl.obolibrary.org/obo/> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX sio:<http://semanticscience.org/resource/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX identifiers:<http://identifiers.org/ensembl/> \n",
            
            "SELECT distinct ?expURI \n",
            "WHERE { \n",            
                "?expURI atlasterms:hasAnalysis ?analysis . \n",   
                "?analysis atlasterms:hasExpressionValue ?value . \n",      
                "?value a ?diffExpType . \n",
                "?diffExpType rdfs:subClassOf atlasterms:DifferentialExpressionRatio . \n",     
                "?value rdfs:label ?expressionValue . \n",     
                "?value atlasterms:pValue ?pvalue . \n",      
                "?value atlasterms:hasFactorValue ?factor . \n",      
                "?value atlasterms:isMeasurementOf ?probe . \n",    
                "?probe atlasterms:dbXref", ensembluri ,". \n",
                
                "}")       
    
    
    
    res <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("an error occured when trying to query for experiment URIs by gene ID ", err)
            })#end tryCatch
    
    return (res$results)  
    
    
}


######
#query for experiments that contain a specific ensembl gene ID and which is reported as diff expressed
######
getExperimentIdsForGeneURI <- function(geneuri, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){    
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX efo:<http://www.ebi.ac.uk/efo/> \n",
            "PREFIX obo:<http://purl.obolibrary.org/obo/> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX sio:<http://semanticscience.org/resource/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX identifiers:<http://identifiers.org/ensembl/> \n",
            
            "SELECT distinct ?experimentid \n",
            "WHERE { \n", 
            "?expUri dcterms:identifier ?experimentid . \n",          
            "?expUri atlasterms:hasAnalysis ?analysis . \n",   
            "?analysis atlasterms:hasExpressionValue ?value . \n",      
            "?value a ?diffExpType . \n",
            "?diffExpType rdfs:subClassOf atlasterms:DifferentialExpressionRatio . \n",     
            "?value rdfs:label ?expressionValue . \n",     
            "?value atlasterms:pValue ?pvalue . \n",      
            "?value atlasterms:hasFactorValue ?factor . \n",      
            "?value atlasterms:isMeasurementOf ?probe . \n",    
            "?probe atlasterms:dbXref", geneuri ,". \n",
            
            "}")       
  
    res <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("an error occured when trying to query for experiment IDs by gene URI ", err)
            })#end tryCatch
    
    return (res$results)  
    
    
}



#########
#get common name for any entity (if available), for example the common gene name 
#requires input parameter of uri of entity in form <http://entityuri> NOTE: including angle brackets
#e.g. cancer for "efo:EFO_0000311"
#########

getLabel <- function(uri, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX efo:<http://www.ebi.ac.uk/efo/> \n",
            "PREFIX obo:<http://purl.obolibrary.org/obo/> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX sio:<http://semanticscience.org/resource/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            
            "SELECT distinct ?label{ \n",    
                uri,"rdfs:label ?label . \n",
                "}")
    
    
    res <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying getLabel query ", err)
            })#end tryCatch
    
    return (res$results)      
}



#########
#query to get pathways associated 
#########
getPathwaysFromGenesAndCondition <- function(condition, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    query <- paste("PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
            "PREFIX sio: <http://semanticscience.org/resource/> \n",
            "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX biopax3: <http://www.biopax.org/release/biopax-level3.owl#> \n",
            
            "SELECT distinct ?pathwayname ?pathway  ?expressionvalue  ?pvalue \n",
            "WHERE { \n",
                "?protein rdf:type biopax3:Protein . \n",
                "?protein biopax3:memberPhysicalEntity \n",
                "[biopax3:entityReference ?dbXref] . \n",
                "?pathway rdf:type biopax3:Pathway . \n",
                "?pathway biopax3:displayName ?pathwayname . \n",
                "?pathway biopax3:pathwayComponent ?reaction . \n",
                "?reaction rdf:type biopax3:BiochemicalReaction . \n",
                "{ \n",
                    "{?reaction ?rel ?protein .} \n",
                    "UNION \n", 
                    "{ \n", 
                        "?reaction ?rel ?complex . \n",
                        "?complex rdf:type biopax3:Complex . \n",
                        "?complex ?comp ?protein . \n",
                        "} \n", 
                    "} \n", 
                "?factor rdf:type ",condition," . \n",
                "?value atlasterms:hasFactorValue ?factor . \n", 
                "?value atlasterms:isMeasurementOf ?probe . \n",
                "?value atlasterms:pValue ?pvalue . \n",
                "?value rdfs:label ?expressionvalue . \n",
                "?probe atlasterms:dbXref ?dbXref . \n",
                "} \n", 
            "ORDER BY ASC (?pvalue) ")

    pathways <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying getPathwaysFromGenesAndCondition query ", err)
            })#end tryCatch
    
    return (pathways$results)
    
}



##########
#draw gene expression levels vs factors for a given Atlas experiment
##########
drawHeatMapForAtlasExperiment <- function(experimentid, tstatsignificance = 5, numoffactorsdiffexpressedacross = 1, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    
    
    experiment <- paste("atlas:",experimentid, sep="")
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
            "PREFIX sio: <http://semanticscience.org/resource/> \n",
            "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#> \n",
            
            
            "SELECT DISTINCT ?genename ?factorLabel ?tStat WHERE { \n",
                experiment,"atlasterms:hasAnalysis ?analysis . \n",     
                "?analysis atlasterms:hasExpressionValue ?value . \n",
                "?value atlasterms:pValue ?pvalue . \n",      
                "?value atlasterms:tStatistic ?tStat . \n",      
                "?value atlasterms:hasFactorValue ?factor . \n",
                "?factor atlasterms:propertyValue ?factorLabel . \n",
                "?value atlasterms:isMeasurementOf ?probe . \n",   
                "?probe atlasterms:dbXref ?dbXref . \n",
                "?dbXref rdfs:label ?genename . \n",
                "} ORDER BY ?genename limit 10000")
    
    message("Executing query... please wait")
    
    d <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying SPARQL query for drawHeatMapForAtlasExperiment function ", err)
            })#end tryCatch
    
    df <- data.frame(Genename=d$results$genename, Factor=factor(d$results$factorLabel), TStat=d$results$tStat, stringsAsFactors=FALSE)
    attach(df)
    genes <- unique(Genename)
    
    # create the matrix for the results and set the row and col names
    values <- matrix(0, nrow=length(genes), ncol=length(unique(Factor)))
    rownames(values) <- genes
    colnames(values) <- unique(Factor)
    
    i<-1
    while (i <= length(Genename)) {
        
        if((df$TStat[i] >= tstatsignificance) || (df$TStat[i] <= (-1*tstatsignificance))){
            
            gn <-df$Genename[i]
            var <-df$Factor[i]
            tstat <-df$TStat[i]
            
            
            rowindex <- match(gn, rownames(values))
            colindex <- match(var, unique(Factor))
            values[rowindex,colindex] <- tstat
        }
        
        i<-i+1
    }
    
    
    
    #fix this
    values<-values[-which(rowSums(values==0) > numoffactorsdiffexpressedacross),]
    #stats::heatmap(values, scale="none", col = cm.colors(256))
    
    par(oma=c(6,2,2,2))
    heatmap(values, scale="none", col = cm.colors(256))
    
    return(values)
}



#########
#get mappings from NCBO for an efo term
#searchuri is an ontology class uri for which efo mappings are to be found
#e.g. <http://purl.bioontology.org/ontology/SNOMEDCT/87163000>  (leukemia in snomed)
#e.g. <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C2985> (diabetes mellitus NCI Thesaurus)
#e.g. <http://purl.bioontology.org/ontology/ICD10CM/J45> (asthma in ICD-10) 
#########
getOntologyMappings <- function(searchuri, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
            "PREFIX sio: <http://semanticscience.org/resource/> \n",
            "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX maps: <http://protege.stanford.edu/ontologies/mappings/mappings.rdfs#> \n",
            
            "SELECT  distinct ?efouri \n",
            "WHERE {  \n",
                "SERVICE <http://sparql.bioontology.org/mappings/sparql?apikey=aa5bfd22-462e-422a-a88c-4055ba36cd1e>{ \n",
                    
                    "?mapping maps:source ",searchuri," ; \n",
                    " maps:target_ontology <http://bioportal.bioontology.org/ontologies/1136> . \n",          
                    "?mapping maps:target ?efouri  \n",        
                    " }\n",
                "}", sep="")
    
    efouri <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying SPARQL query for getOntologyMappings function ", err)
            })#end tryCatch
    
    return(efouri$results)   
}



#######
#get genes for a given pubmedid, if the experiment is in Atlas
#e.g. "11027337"
######
getGeneListFromPubmedid <- function(searchid, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
            
    searchuri<- paste("<http://identifiers.org/pubmed/",searchid,">", sep="")
            
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX efo:<http://www.ebi.ac.uk/efo/> \n",
            "PREFIX obo:<http://purl.obolibrary.org/obo/> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX sio:<http://semanticscience.org/resource/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            
            "SELECT distinct ?experimenturi \n",       
            "WHERE {\n",           
            "?experimenturi rdf:type atlasterms:Experiment . \n",  
            "?experimenturi atlasterms:pubmedid", searchuri ," . \n",        
            "}")
    
    message("Performing query please wait...")

    
    res <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying SPARQL query for getGeneListFromPubmedid function ", err)
            })#end tryCatch
    
    message("Getting genes for experiment ", res$results[1])

    genelist <- getGenesForExperimentURI(res$results[1], endpoint)    
    
    return(genelist)      
 
}


##########
#get Reactome pathway for Ensembl gene id
##########
getGenesForPathwayURI <- function(pathwayuri, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    
    
    query <- paste("PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
            "PREFIX sio: <http://semanticscience.org/resource/> \n",
            "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#> \n",
            "PREFIX biopax3:<http://www.biopax.org/release/biopax-level3.owl#> \n",
            
            
            "SELECT distinct ?geneuri \n",
            "WHERE { \n",
                
                
                # query for pathways by those protein targets
                "SERVICE <http://www.ebi.ac.uk/rdf/services/reactome/sparql> { \n",                            
                    pathwayuri ,"rdf:type biopax3:Pathway .  \n",
                    pathwayuri ,"biopax3:pathwayComponent ?reaction . \n",
                    "?reaction rdf:type biopax3:BiochemicalReaction . \n",
                    "{ \n",         
                        "{?reaction ?rel ?protein .} \n",  
                        "UNION  \n",
                        "{ \n", 
                            "?reaction  ?rel  ?complex . \n",
                            "?complex rdf:type biopax3:Complex . \n", 
                            "?complex ?comp ?protein . \n",
                            "}} \n", 
                    "?protein rdf:type biopax3:Protein . \n",
                    "?protein biopax3:memberPhysicalEntity \n",
                    "[biopax3:entityReference ?dbXref ] . \n",
                    "} \n",   
                
                # get Atlas experiment plus experimental factor where protein is expressed

                "?probe atlasterms:dbXref ?dbXref . \n",
                "?probe atlasterms:dbXref ?geneuri . \n",
                "?geneuri rdf:type atlasterms:EnsemblDatabaseReference . \n",
              
                               
                "} \n"
            )  #endpaste   

    
    genes <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying SPARQL query for getGenesForPathwayURI function ", err)
            })#end tryCatch
  
    return(genes$results)   
    
}




##########
#get Reactome pathway for Ensembl gene id
##########
getPathwayForGeneId <- function(geneid,  endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
   
    #get gene uri
    geneuri <- getGeneUriFromEnsemblId(geneid)
    
    message(geneuri)
    
    if(length(geneuri) != 0){
        
        
        query <- paste("PREFIX atlas_r: <http://atlasrdfrpackage> \n",
                "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
                "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
                "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
                "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
                "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
                "PREFIX sio: <http://semanticscience.org/resource/> \n",
                "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
                "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
                "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
                "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#> \n",
                "PREFIX biopax3:<http://www.biopax.org/release/biopax-level3.owl#> \n",
                
                
                "SELECT distinct ?pathway(str(?pathwayname) as ?pathname) \n",
                "WHERE { \n",
                    
                    # get Atlas experiment plus experimental factor where protein is expressed
                    "?probe atlasterms:dbXref ", geneuri ," . \n",
                    "?probe atlasterms:dbXref ?dbXref . \n",             
                    
                    # query for pathways by those protein targets
                    "SERVICE <http://www.ebi.ac.uk/rdf/services/reactome/sparql> { \n",                            
                        "?pathway rdf:type biopax3:Pathway .  \n",
                        "?pathway ?p ?pathwayname . \n",
                        "?pathway biopax3:pathwayComponent ?reaction . \n",
                        "?reaction rdf:type biopax3:BiochemicalReaction . \n",
                        "{ \n",         
                            "{?reaction ?rel ?protein .} \n",  
                            "UNION  \n",
                            "{ \n", 
                                "?reaction  ?rel  ?complex . \n",
                                "?complex rdf:type biopax3:Complex . \n", 
                                "?complex ?comp ?protein . \n",
                                "}} \n", 
                        "?protein rdf:type biopax3:Protein . \n",
                        "?protein biopax3:memberPhysicalEntity \n",
                        "[biopax3:entityReference ?dbXref ] . \n",
                        "} \n",  
                    
                    "filter regex(str(?p), \"displayName\") . \n", 
                    
                    "} \n"
                )  #endpaste   
        
        
        
        pathways <- tryCatch({
                    SPARQL(url=endpoint,query)
                },
                error = function(err){
                    message("An error occured when trying SPARQL query for getPathwayForGeneId function ", err)
                })#end tryCatch
        
        return(pathways$results)
        
    }
    else{
        message("Could not find gene ", geneid)
        
    }
    
    
}




##########
#getRankedPathwaysForGeneList - function to get pathways which are associated to a gene list, as ranked by pathways with most genes associated
##########

getRankedPathwaysForGeneIds <- function(genelist, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){

    
    rankedpathways <- list()
    
    #loop through list
    for (i in 1:length(genelist)){
        
        #get gene uri
        geneuri <- getGeneUriFromEnsemblId(genelist[i])
        
        if(length(geneuri) != 0){
            
            message("querying ", geneuri)
            
            
            query <- paste("PREFIX atlas_r: <http://atlasrdfrpackage> \n",
                    "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
                    "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
                    "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
                    "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
                    "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
                    "PREFIX sio: <http://semanticscience.org/resource/> \n",
                    "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
                    "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
                    "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
                    "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#> \n",
                    "PREFIX biopax3:<http://www.biopax.org/release/biopax-level3.owl#> \n",
                    
                    
                    "SELECT distinct ?pathway(str(?pathwayname) as ?pathname) \n",
                    "WHERE { \n",
                        
                        # get Atlas experiment plus experimental factor where protein is expressed
                        "?probe atlasterms:dbXref ", geneuri ," . \n",
                        "?probe atlasterms:dbXref ?dbXref . \n",             
                        
                        # query for pathways by those protein targets
                        "SERVICE <http://www.ebi.ac.uk/rdf/services/reactome/sparql> { \n",                            
                            "?pathway rdf:type biopax3:Pathway .  \n",
                            "?pathway ?p ?pathwayname . \n",
                            "?pathway biopax3:pathwayComponent ?reaction . \n",
                            "?reaction rdf:type biopax3:BiochemicalReaction . \n",
                            "{ \n",         
                                "{?reaction ?rel ?protein .} \n",  
                                "UNION  \n",
                                "{ \n", 
                                    "?reaction  ?rel  ?complex . \n",
                                    "?complex rdf:type biopax3:Complex . \n", 
                                    "?complex ?comp ?protein . \n",
                                    "}} \n", 
                            "?protein rdf:type biopax3:Protein . \n",
                            "?protein biopax3:memberPhysicalEntity \n",
                            "[biopax3:entityReference ?dbXref ] . \n",
                            "} \n",  
                        
                        "filter regex(str(?p), \"displayName\") . \n", 
                        
                        "} \n"
                    )  #endpaste   

           
            pathways <- tryCatch({
                        SPARQL(url=endpoint,query)
                    },
                    error = function(err){
                        message("An error occured when trying SPARQL query for getRankedPathwaysForGeneIds function ", err)
                    })#end tryCatch

            
            result <- pathways$results
            
            foundpathway = FALSE
            if(length(result) != 0){
                #for each pathway returned
                for(j in 1:nrow(result)){
                    
                    #if this isn't first time through the list, otherwise the list will be empty
                    if(length(rankedpathways) !=0 ){
                        
                        #see if the pathway is already in list just add genes to it
                        for(k in 1:length(rankedpathways)){
                            
                            if(result[j,1] == rankedpathways[[k]]@pathwayuri){
                                
                                #add genes to vector  
                                vectorgenes <- rankedpathways[[k]]@genes
                                vectorgenes <- c(vectorgenes, geneuri)
                                rankedpathways[[k]]@genes <- vectorgenes
                                rankedpathways[[k]]@numgenes <- length(vectorgenes)
                                
                                #since it exists break out of loop
                                foundpathway = TRUE
                                break
                            }                                              
                        }
                        
                    }    
                    
                    #if the pathway hasn't been found then add a new one
                    if(foundpathway == FALSE){
                        
                        singlepathway <- new("pathwayresult")
                        
                        singlepathway@pathwayuri <- result[j,1]  #pathway uri
                        singlepathway@label <- result[j,2]  #pathway name
                        singlepathway@genes <- c(geneuri)  #gene for this pathway
                        singlepathway@numgenes <- 1 # set gene count to 1
                        rankedpathways <- c(rankedpathways, singlepathway)
                    }
                    
                    
                }
                
            }#end if 
        }
        
        
    }#end for loop sparql query    
    
    if(length(rankedpathways) > 1){
        #mint a blank matrix
        vec <- c(0,0)
        matrixcounts <- matrix()
        matrixcounts <- rbind(vec)
        
        #order results by pathway with most genes
        for(i in 1:length(rankedpathways)){
            
            #extract out the gene counts
            vec <- c(i, rankedpathways[[i]]@numgenes)
            matrixcounts <- rbind(matrixcounts, vec)
            
        }
        
        #get rid of blank row at top
        matrixcounts <- matrixcounts[-1,]
        #order by gene counts
        matrixcounts <- matrixcounts[order(matrixcounts[,2], decreasing=TRUE),]
        
        sortedpathways <- rankedpathways
        #for each matrix count, reorder the pathway results with largest first
        for(i in 1:nrow(matrixcounts)){
            
            sortedpathways[[i]] <- rankedpathways[[matrixcounts[i,1]]]             
        }        
    }
    else {
        sortedpathways <- rankedpathways    
    }
    
    return(sortedpathways)    
}



###########
#get efo URI based on label - exact match only
###########
searchForEFOTerms <- function(label, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
            "PREFIX sio: <http://semanticscience.org/resource/> \n",
            "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#> \n",
            
            "SELECT distinct ?efouri ?label WHERE { \n",
            "?efouri rdfs:subClassOf* efo:EFO_0000001 . \n",
            "?efouri rdfs:label ?label . \n", 
            "FILTER regex(str(?label), \"",label,"\", \"i\"). \n", 
            "}", sep="")
    
    
    uris <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying SPARQL query for searchForEFOTerms function ", err)
            })#end tryCatch
    
    
    uris <- SPARQL(url=endpoint, query)
    return (uris$results)
    
    
}


###################
#functions to perform enrichment analysis
###################


###########
#get gene uris given ensembl gene name and a taxon
###########
getGeneUriFromName <- function(genename, taxon, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
            "PREFIX sio: <http://semanticscience.org/resource/> \n",
            "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#> \n",
            
            "SELECT distinct ?geneuri WHERE { \n",
                "?geneuri rdf:type atlasterms:EnsemblDatabaseReference . \n",
                "?geneuri atlasterms:taxon ", taxon , ". \n",
                "?geneuri rdfs:label ?label. \n", 
                "FILTER regex(str(?label), \"^",genename,"$\", \"i\"). \n", 
                "}", sep="")
    
    uris <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying SPARQL query for getGeneUriFromName function ", err)
            })#end tryCatch
    
    return (uris$results)
    
}

###########
#get efo URI based on label - exact match only
###########
getExFactorURIFromLabel <- function(label, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    query <- paste("PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
            "PREFIX sio: <http://semanticscience.org/resource/> \n",
            "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#> \n",
            
            "SELECT distinct ?efouri WHERE { \n",
                "?efouri rdfs:subClassOf* efo:EFO_0000001 . \n",
                "?efouri rdfs:label ?label . \n", 
                "FILTER regex(str(?label), \"^",label,"$\", \"i\"). \n", 
                "}", sep="")

    
    uris <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying SPARQL query for getExFactorURIFromLabel function ", err)
            })#end tryCatch
    return (uris$results)
    
    
    
}


###########
#get efo URI based on label - exact match only
###########
searchForEFOTerms <- function(label, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
            "PREFIX sio: <http://semanticscience.org/resource/> \n",
            "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#> \n",
            
            "SELECT distinct ?efouri ?label WHERE { \n",
            "?efouri rdfs:subClassOf* efo:EFO_0000001 . \n",
            "?efouri rdfs:label ?label . \n", 
                "FILTER regex(str(?label), \"",label,"\", \"i\"). \n", 
            "}", sep="")

    
    uris <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying SPARQL query for searchForEFOTerms function ", err)
            })#end tryCatch
    return (uris$results)
    
    
}




###########
#get gene uris given ensembl gene id 
###########
getGeneUriFromEnsemblId <- function(id, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
            "PREFIX sio: <http://semanticscience.org/resource/> \n",
            "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#> \n",
            
            "SELECT distinct ?geneuri WHERE { \n",
                "?geneuri rdf:type atlasterms:EnsemblDatabaseReference . \n",
                "?geneuri dcterms:identifier ?geneid. \n", 
                "FILTER regex(str(?geneid), \"",id,"\", \"i\"). \n", 
                "}", sep="")
    
    uris <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying SPARQL query for getGeneUriFromEnsemblId function ", err)
            })#end tryCatch
    return (uris$results)
    
}


###########
#get pathway uri for given pathway name
###########
getPathwayUriFromName <- function(name, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    
    
    query <- paste( "PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
            "PREFIX sio: <http://semanticscience.org/resource/> \n",
            "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#> \n",
            "PREFIX biopax3:<http://www.biopax.org/release/biopax-level3.owl#> \n",
            
            "SELECT distinct ?pathwayuri ?label WHERE { \n",
                "?pathwayuri rdf:type biopax3:Pathway . \n",
                "?pathwayuri biopax3:displayName ?label . \n",
                "FILTER regex(str(?label), \"",name,"\", \"i\"). \n", 
                "}", sep="")
    
    uris <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying SPARQL query for getPathwayUriFromName function ", err)
            })#end tryCatch
    return (uris$results)
    
}


#######
#convenience function to get gene counts for the gene list in question
#genelist should be a vector of gene uris which are being analysed for enrichment
#######
calculateCountsForGeneLists <- function(genelist, genelist_bg, genecounts){

    
    #creat a hash to store the genelist's factors in (hash of factorbackground classes)
    genelistfactors <- hash()
    
    #preprocess the genelist to check for multiple uris for a single gene - this needs to be flattened
    processedgenelist <- vector()
    for(i in 1:length(genelist)){   
        if(length(genelist[[i]]) > 1){
            for(j in 1:length(genelist[[i]])){
                processedgenelist <- c(processedgenelist, genelist[[i]][j])
            }
        }
        else{
            processedgenelist <- c(processedgenelist, genelist[[i]])
        }
        
    }
    
    #replace the gene list with the normalised flat gene list
    genelist <- processedgenelist
    #get total genes for this gene list 
    totalnumbergenes <- length(genelist)
    
    #for each gene in gene list find the ex factors terms associated  
    for(i in 1:length(genelist)){    
        
        geneobject <- genelist_bg[[genelist[[i]]]]
        
        #if the gene has been found in the background
        if(!is.null(geneobject)){
            message("Found gene ", genelist[[i]] )
            #get the ex factors for this gene
            singlegeneexfactors <- geneobject$exfactoruris
            species <- geneobject$species
            
            #for each ex factor for this gene 
            for(j in 1:length(singlegeneexfactors)){
                #check to see if the factor exists already
                if(has.key(singlegeneexfactors[j], genelistfactors) == TRUE){      
                    #the ex factor is in the hash set already so merge gene into slot
                    factorobject <- genelistfactors[[singlegeneexfactors[j]]]                  
                    factorobject$mergegeneuris(genelist[[i]])
                    numgenesex <- length(factorobject$getgeneuris())
                    factorobject$setnumgenesexpressed(as.integer(numgenesex))
                    factorobject$setnumgenesnotexpressed(as.integer(totalnumbergenes-numgenesex))
                    factorobject$species <- species                    
                    
                    #add back to the hash set
                    .set(genelistfactors, keys=singlegeneexfactors[j], values=factorobject)
                }
                #otherwise mint new object and add to hash
                else{
                    factorobject <- factorbackground$new()
                    factorobject$seturi(singlegeneexfactors[j])
                    factorobject$mergegeneuris(genelist[[i]])
                    numgenesex <- length(factorobject$getgeneuris())
                    factorobject$setnumgenesexpressed(as.integer(numgenesex))
                    factorobject$setnumgenesnotexpressed(as.integer(totalnumbergenes-numgenesex))
                    factorobject$species <- species
                    
                    #add new object to the hash set
                    .set(genelistfactors, keys=singlegeneexfactors[j], values=factorobject)
                    
                }
            }#end for             
        }
        else{
            message("Did not find gene ", genelist[[i]] , ": ignoring")          
        }
    }
    
    return(genelistfactors)
    
}



###############
#function to do enrichment using Fishers exact test based on a set of genelistfactors (factors for genes) and bg counts
#the input requires the gene list to be in the form of identifiers.org ensembl uris to turn common names into
#uris use the function getEnsemblUrisFromNames
###############
doFishersEnrichment <- function(genelist, genelist_bg, genecounts){
    
    #calc counts for given gene list
    genelistfactors <- calculateCountsForGeneLists(genelist, genelist_bg, genecounts)
    
    
    #specify class to store enrichemt results 
    enrichmentresult <- setClass("enrichmentresult",           
            representation( factoruri="character", 
                    label="character", 
                    p.value="numeric", 
                    estimate="numeric",    
                    alternative="character", 
                    null.value="numeric", 
                    method="character",
                    enrichedgenes="vector"))
    
    fisherresults <- list()
    
    #do test for each factor
    genes <- keys(genelistfactors)
    for (i in 1:length(genes)){
        
        genelistobject <- genelistfactors[[genes[i]]]
        bggeneobject <- genecounts[[genes[i]]]
        
        if(!is.null(genelistobject) && !is.null(bggeneobject)){
            
            #gather stats for fisher test
            genelistannotated <- genelistobject$numgenesexpressed
            genelistnotannotated <- genelistobject$numgenesnotexpressed
            bgannotated <- bggeneobject$numgenesexpressed
            bgnotannotated <- bggeneobject$numgenesnotexpressed
            
            input <- matrix(c(genelistannotated, genelistnotannotated, bgannotated, bgnotannotated), nrow = 2, dimnames =
                    list(c("Annotated", "Not Annotated"),
                            c("Genelist", "Backgound")))
            
            #do fisher's exact test
            result <- fisher.test(input)
            
            #store results
            enrichresult <- new("enrichmentresult")
            
            enrichresult@label <- bggeneobject$label
            enrichresult@factoruri <- genelistobject$uri
            enrichresult@enrichedgenes <- genelistobject$geneuris
            enrichresult@p.value <- result$p.value            
            enrichresult@estimate <- result$estimate
            enrichresult@alternative <- result$alternative
            enrichresult@null.value <- result$null.value
            enrichresult@method <- result$method
            
            fisherresults <- c(fisherresults, enrichresult)                        
            
        }
        
    }
    message("enrichment complete")
    return(fisherresults)
    
}




###############
#function to do enrichment analysis using Fisher's exact test based on a gene list specified by common gene name
#input requires a vector of the common gene names and a taxon id e.g. obo:NCBITaxon_9606 for human (note genes from multiple species not allowed) 
###############
doFishersEnrichmentForGeneNames <- function(genenames, taxon, genelist_bg, genecounts, endpoint="http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    geneuris <- vector()
    
    for(i in 1:length(genenames)){
        
        geneuris <- c(geneuris, getGeneUriFromName(genenames[[i]], taxon, endpoint))
    }
    
    #if none of these gens are found then do not go further
    if(length(geneuris)==0){
        message("None of the genelist was found in the Atlas data. Halting enrichment.")
    }    
    else{
        results <- doFishersEnrichment(geneuris, genelist_bg, genecounts)
        return(results)
    }
}




###############
#function to do enrichment analysis using Fisher's exact test based on a gene list specified by common gene name
#input requires a vector of the common gene names
###############
doFishersEnrichmentForEnsemblIds <- function(geneids, genelist_bg, genecounts){
    
    geneuris <- vector()
    
    for(i in 1:length(geneids)){
        
        geneuris <- c(geneuris, getGeneUriFromEnsemblId(geneids[[i]]))
    }
    
    #if none of these gens are found then do not go further
    if(length(geneuris)==0){
        message("None of the genelist was found in the Atlas data. Halting enrichment.")
    }   
    else{
        results <- doFishersEnrichment(geneuris, genelist_bg, genecounts)
        return(results)
    }
}



##############
#function to filter set of enriched factors given a specific p-value cutoff, default is 0.05
#results are then vizualised as a bar plot and also returned as a vector of pvalues to factor
##############
vizPvalues <- function(resultset, cutoff = "0.05"){
    
    
    results <- vector()
    names <- vector()
    
    for (i in 1:length(resultset)){
        
        pvalue <- resultset[[i]]@p.value
        #if the p value is below the user cut off
        if (round(pvalue, digits=5) <= cutoff){
            #message(resultset[[i]]@label, "   ", pvalue)     
            
            results <- c(results, pvalue)
            names <- c(names, resultset[[i]]@label)
            results <- setNames(results, names)    
            
        }  
        
    }
    
    results <- sort(results)
    
    par(mar=c(3,12,1,1))
    barplot(height=results, names.arg=names(results), horiz = TRUE, las=1, xlab = "p-value", 
            col = hcl(seq(0, 200, length = length(results))), cex.names=0.8)
    title(main = list(paste("Enrichment with p-vlaue cut off ",cutoff ), font = 3))
    
    return(results)
}

###########
#function to filter out subclasses of a given class from the enrichment results
#e.g. filterparentclass="obo:CHEBI_37577" would filter out all chemical compounds from the results
###########
excludeSubclasses <- function(filterparentclass, resultset, endpoint="http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    #get the subclass uris for the class to filter out
    query <- paste("PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
            "PREFIX sio: <http://semanticscience.org/resource/> \n",
            "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#> \n",
            
            "SELECT distinct ?classuri WHERE { \n",           
            "?classuri rdfs:subClassOf*", filterparentclass,". \n",    
            "}")
    
    uris <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying SPARQL query for excludeSubclasses function ", err)
            })#end tryCatch
    
    sizeunfilitered <- length(resultset)     
    
    #create vector of uris in the result set
    resulturis <- vector()
    for(i in 1:length(resultset)){
        resulturis <- c(resulturis, resultset[[i]]@factoruri)
    }
    
    
    #loop through each result and remove it if it exists in the results
    for(i in 1:nrow(uris$results)){
        
        #if this is not null then there is a match
        if(uris$results[i,] %in% resulturis){
            
            #get element number but from filtered list as this will be different because it's being resized
            removelement <- match(uris$results[i,], resulturis)
            #remove it from the list
            resultset <- resultset[-(removelement)]
            #and from vector
            resulturis <- resulturis[-(removelement)]     
        }       
    }  
    message("Removed " ,(sizeunfilitered - length(resultset))," factors from result set")    
    
    return (resultset)    
}




###########
#function to include only subclasses of a given class from the enrichment results, removing all the rest
#e.g. includeparentclass="obo:CHEBI_37577" would filter out all chemical compounds from the results
###########
includeOnlySubclasses <- function(includeparentclass, resultset, endpoint="http://www.ebi.ac.uk/rdf/services/atlas/sparql"){
    
    #create new filtered results
    filteredresults <- vector()
    
    #get the subclass uris for the class to filter out
    query <- paste("PREFIX atlas_r: <http://atlasrdfrpackage> \n",
            "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> \n",
            "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> \n",
            "PREFIX owl: <http://www.w3.org/2002/07/owl#> \n",
            "PREFIX dcterms: <http://purl.org/dc/terms/> \n",
            "PREFIX obo: <http://purl.obolibrary.org/obo/> \n",
            "PREFIX sio: <http://semanticscience.org/resource/> \n",
            "PREFIX efo: <http://www.ebi.ac.uk/efo/> \n",
            "PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/> \n",
            "PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/> \n",
            "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#> \n",
            
            "SELECT distinct ?classuri WHERE { \n",           
                "?classuri rdfs:subClassOf*", includeparentclass,". \n",    
                "}")
    
    
    uris <- tryCatch({
                SPARQL(url=endpoint,query)
            },
            error = function(err){
                message("An error occured when trying SPARQL query for includeOnlySubclasses function ", err)
            })#end tryCatch
    
    #if subclasses were successfully retrieved from the sparql query
    if(length(uris$results > 0)){
        
        #create vector of uris in the result set
        resulturis <- vector()
        for(i in 1:length(resultset)){
            resulturis <- c(resulturis, resultset[[i]]@factoruri)
        }
        
        #loop through each result and remove it if it exists in the results
        for(i in 1:length(uris$results)){
            
            #if this is not null then there is a match
            if(uris$results[i] %in% resulturis){
                
                #get element number but from filtered list as this will be different because it's being resized
                includeelement <- match(uris$results[i], resulturis)
                #remove it from the list
                filteredresults <- c(filteredresults,resultset[includeelement])         
            }        
        }
        message("Removed " ,(length(resulturis) - length(filteredresults))," factors from result set")    
        
        return (filteredresults)  
        
    }
    #this may be because the sparql endpoint is down or the subclass given does not have any child classes
    else{
        message("Error in retrieving subclasses of class ", includeparentclass)
    }
}



###########
#function to get uri of species based on name. Used for functions where species can be specified 
#works for human, rat, mouse, arabidopsis and drosophila
###########
getTaxonURI <- function(taxonName){
    
    taxonName <- tolower(taxonName)
    
    if(taxonName == "human" || taxonName == "homo sapiens"){
        message("human")
        uri <- "obo:NCBITaxon_9606"
        
    }
    else if(taxonName == "mouse" || taxonName == "mus musculus"){
        message("mouse")
        uri <- "obo:NCBITaxon_10090"
    }
    else if(taxonName == "arabidopsis thaliana" || taxonName == "arabidopsis"){
        message("arabidopsis")
        uri <- "obo:NCBITaxon_3702"
    }   
    else if(taxonName == "rat" || taxonName == "rattus norvegicus"){
        message("rattus norvegicus")
        uri <- "obo:NCBITaxon_10116"
    }
    else if(taxonName == "fly" || taxonName == "drosophila" || taxonName == "drosophila melanogaster"){
        message("drosophila")
        uri <- "obo:NCBITaxon_7227"
    }
    else{
        message("Could not identify species")
    }
    
    return(uri)
    
}


###########
#order enrichment results set by p value
###########
orderEnrichmentResults <- function(resultset){
    
    sortedresultset <- list()
    
    resultpvalues <- data.frame(originalposition=character(), pvalue=character())
    
    #extract the pvlues from the result set
    for(i in 1:length(resultset)){
        newrow <- c(i, resultset[[i]]@p.value)
        resultpvalues <- rbind(resultpvalues,newrow)
    }

    sortedpvalues <- resultpvalues[ order(resultpvalues[,2]),]  
    
    for(i in 1:nrow(sortedpvalues)){       
        sortedresultset <- c(sortedresultset, resultset[[sortedpvalues[i,1]]])        
    }
   
    return(sortedresultset)   
    
}


