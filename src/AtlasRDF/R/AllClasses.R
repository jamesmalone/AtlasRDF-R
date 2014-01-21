#All classes used in the package


#class to store information about a pathway and associated genes
setClass("pathwayresult",           
        representation( pathwayuri="character", 
                label="character",    
                numgenes="numeric", 
                genes="vector"))



#class used to store gene references, including the URI, name, ID, species and experimental factors associated with the gene
generef <- setRefClass("generef",
        fields = list( geneuri = "character",
                genelabel = "character",
                geneensemblid = "character",
                species = "character",
                exfactoruris="vector"),
        
        methods = list(
                setgeneuri = function(value) {
                    geneuri <<- value
                    invisible(value)
                },                
                getgeneuri = function() {
                    return(geneuri)
                },                
                setgenelabel = function(value) {
                    genelabel <<- value
                    invisible(value)
                },
                getgenelabel = function() {
                    return(genelabel)
                },
                setensemblid = function(value) {
                    geneuri <<- value
                    invisible(value)
                },                
                getensemblid = function() {
                    return(geneuri)
                },  
                setspecies = function(value) {
                    species <<- value
                    invisible(value)
                },                
                getspecies = function() {
                    return(species)
                },  
                getexfactoruris = function() {
                    return(exfactoruris)
                }, 
                show = function() {
                    message("An object of class generef used to store gene reference objects")
                    message("Field values")
                    message("geneuri: ", geneuri) 
                    message("genelabel: ", genelabel)
                    message("geneensemblid: ", geneensemblid)
                    message("species: ", species)
                    message("exfactoruris (URI identifiers for experimental factors associated with this gene): ", exfactoruris)
                   
                },   
                #function to merge a passed parameter with the existing set of genes stored in this object
                mergeexfactoruris = function(value){
                    if(!is.null(value)){
                        exfactoruris <<-unique(c(exfactoruris, value))
                    }                      
                }        
                ))



#the factor background class used in enrichment to represent an experimental factor and any associated genes for that 
#factor. also stores info on super and subclasses of this factor
factorbackground <- setRefClass("factorbackground",
        fields = list( uri = "character",
                label = "character",
                species = "character",
                geneuris="vector",
                numgenesexpressed="integer", 
                numgenesnotexpressed="integer", 
                subclasses="vector", 
                superclasses="vector"),
        
        methods = list(
                seturi = function(value) {
                    uri <<- value
                    invisible(value)
                },
                
                geturi = function() {
                    return(uri)
                },
                
                setlabel = function(value) {
                    label <<- value
                    invisible(value)
                },
                
                getlabel = function() {
                    return(label)
                },
                
                setspecies = function(value) {
                    species <<- value
                    invisible(value)
                },
                
                getspecies = function() {
                    return(species)
                },
                
                setgeneuris = function(value) {
                    geneuris <<- value
                    invisible(value)
                },
                
                getgeneuris = function() {
                    return(geneuris)
                },
                
                setnumgenesexpressed = function(value) {
                    numgenesexpressed <<- value
                    invisible(value)
                },
                
                getnumgenesexpressed = function() {
                    return(numgenesexpressed)
                },
                
                setnumgenesnotexpressed = function(value) {
                    numgenesnotexpressed <<- value
                    invisible(value)
                },
                
                getnumgenesnotexpressed = function() {
                    return(numgenesnotexpressed)
                },
                
                setsubclasses = function(value) {
                    subclasses <<- value
                    invisible(value)
                },
                
                getsubclasses = function() {
                    return(subclasses)
                },
                
                setsuperclasses = function(value) {
                    superclasses <<- value
                    invisible(value)
                },
                
                getsuperclasses = function() {
                    return(superclasses)
                },
                show = function() {
                    message("An object of class factorbackground used to store an experimental factor and any gene associated with them")
                    message("Field values")
                    message("uri (identifier for the factor): ", uri) 
                    message("label (human readable label for the factor): ", label) 
                    message("species (the species of the genes associated with the factor): ", species) 
                    message("geneuris (the URI identifiers of the genes associated with the factor): ", geneuris) 
                    message("numgenesexpressed (the total number of genes differetially expressed and associated with the factor): ", numgenesexpressed)
                    message("numgenesnotexpressed (the total number of genes not differetially expressed and associated with the factor): ", numgenesnotexpressed)  
                    message("subclasses (ontology subclasses for this factor): ", subclasses) 
                    message("superclasses (ontology superclasses for this factor): ", superclasses) 
                }, 
                
                #function to merge a passed parameter with the existing set of genes stored in this object
                mergegeneuris = function(value){
                    if(!is.null(value)){
                        geneuris <<-unique(c(geneuris, value))
                    } 
                }
                
                ))



#class to store enrichemt results, including the uri of the factor plus information on the statistics for enrichment
#such as p.value and the genes enriched for this factor
setClass("enrichmentresult",           
        representation( factoruri="character", 
                label="character", 
                p.value="numeric", 
                estimate="numeric",    
                alternative="character", 
                null.value="numeric", 
                method="character",
                enrichedgenes="vector"))




