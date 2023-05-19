#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyalert)
library(shinylogs)
library(seqinr)
library(uuid)
# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    
    track_usage(storage_mode = store_json(path = "logs/"))
    server_config = read.table("/home/stotoshka/server_config.txt",sep="=",row.names = 1)
    uuid = UUIDgenerate(output = "string")
    input.dst = file.path(server_config["userdata_path","V2"],paste0(uuid,".txt"))
    
    write.input.to.file = function(sequences){
        fileConn<-file(input.dst)
        writeLines(sequences, fileConn)
        close(fileConn)
    }
    
    check.input = function(sequences){
        all(grepl("^[ACDEFGHIKLMNPQRSTVWY]+$",sequences))
    }
    
    observeEvent(input$apply,{
        tryCatch({
            if (input$input.text != ""){
                sequences = input$input.text 
                write.input.to.file(sequences)
            }
            #Проверить вход
            input.seqs = readLines(input.dst)
            if (check.input(input.seqs)){
                shinyalert("Debug",
                           "Good input",
                           type = "info"
                )
            }
            else{
                stop("Incorrect input. Sequences must be written in one-letter code.")
            }
        },
        error = function(e){
            shinyalert("Error",
                       paste("Не могу считать данные, потому что",e),
                       type = "error"
            )
        })
    })
    
    observeEvent(input$reset,{
        session$reload()
        shinyalert("Debug",
                   as.character(file.exists(input.dst)),
                   type = "info"
        )
    })
    
    session$onSessionEnded(function() {
        file.remove(input.dst)
        print(paste(uuid, 'the session has ended'))
    })
})
