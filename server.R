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
library(stringr)
library(dplyr)
library(promises)
library(ssh)
library(future)
library(uuid)
Sys.setlocale(category = "LC_ALL", locale = "ru_RU.UTF-8")
# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    
    server_config = read.table("/home/stotoshka/server_config.txt",sep="=",row.names = 1)
    pass_server = server_config["pass_server","V2"]
    pass_passwd = server_config["pass_server_password","V2"]
    uuid = str_replace_all(UUIDgenerate(output = "string"),"-","_")
    input.dst = file.path(server_config["userdata_path","V2"],paste0(uuid,".csv"))
    run.logs = file.path(server_config["userdata_path","V2"],"logs",paste0(uuid,".log"))
    track_usage(storage_mode = store_json(path = file.path(server_config["userdata_path","V2"],"shiny")))
    
    
    final.result = reactiveValues(
        proteasome = NULL,
        TAP = NULL,
        HLA = NULL
    )


    write.input.to.file = function(sequences){
        fileConn<-file(input.dst)
        writeLines(c("proteins",sequences), fileConn)
        close(fileConn)
    }
    
    check.input = function(sequences){
        all(grepl("^[ACDEFGHIKLMNPQRSTVWY]+$",sequences))
    }
    
    get_status <- function(){
        scan(run.logs, what = "character",sep="\n")
    }
    
    set_status <- function(msg){
        write(paste(msg, format(Sys.time(), "%a %b %d %X %Y")), run.logs,append = T)
    }
    
    launch.cmd.with.waiting = function(ssh, command,success_status){
        out = ssh_exec_internal(ssh, command)
        if (out$status == 0){
            set_status(success_status)
        }
        else{
            set_status("Error")
            set_status(rawToChar(out$stdout))
            set_status(rawToChar(out$stderr))
            stop("Error. See log.")
        }
    }
    
    observeEvent(input$apply,{
        tryCatch({
            if (input$input.text != ""){
                sequences = input$input.text 
                write.input.to.file(sequences)
            }
            #Проверить вход
            input.seqs = read.csv2(input.dst)$proteins
            if (check.input(input.seqs)){
                disable("apply")
                disable("reset")
                ssh.session = ssh_connect(pass_server, passwd = pass_passwd)

                result <- future({
                    set_status("Запуск")
                    workdir = paste0(server_config["pass_root","V2"],"\\", uuid)
                    ssh_exec_wait(ssh.session, paste("mkdir",workdir))
                    scp_upload(ssh.session, input.dst, to = paste0(uuid,".csv"))#WARNING Cannot directly move to needed folder
                    ssh_exec_wait(ssh.session, paste("move",paste0("'",uuid,".csv'"), workdir))
                    ssh_exec_wait(ssh.session, paste("cd",workdir))
                    launch.cmd.with.waiting(ssh.session,command = paste0("python ",server_config["python_scripts","V2"],"\\renameInput.py ",workdir),
                                            success_status = "Входные последовательности получены.")#Убрать эти кавычки странные
                    remote.input.path = paste0(workdir,"\\", uuid, ".csv")
                    path_to_config = paste0(workdir,"\\", uuid, ".json")
                    launch.cmd.with.waiting(ssh.session,command = paste0("python ",server_config["python_scripts","V2"],"\\generateConfigForConverter.py -i ",
                                                                         remote.input.path," -o ", workdir," -c proteins -u -t 5 ",path_to_config),
                                            success_status = "Есть конфигурация для конвертера.")
                    launch.cmd.with.waiting(ssh.session,command = paste0("python ",server_config["python_scripts","V2"],"\\SeqToSDF.py ",path_to_config),
                                            success_status = "Сконвертировано в SDF.")
                   
                    final.result$proteasome = data.frame(matrix(1:10,nrow = 5))
                    final.result$TAP = data.frame(matrix(1:10,nrow = 5))
                    final.result$HLA = data.frame(matrix(1:10,nrow = 5))
                })
                
                # Catch inturrupt (or any other error) and notify user
                result <- catch(result,
                                function(e){
                                    final.result(NULL)
                                    set_status(paste("Ошибка",e$message))
                                })
                
                # After the promise has been evaluated set nclicks to 0 to allow for anlother Run
                result <- finally(result,
                                  function(){
                                      tryCatch({
                                          file.remove(input.dst)
                                          #launch.cmd.with.waiting(ssh.session,command = paste0("python ",server_config["python_scripts","V2"],"\\deleteRemote.py ",uuid),
                                          #                    success_status = "Очищено.")
                                          ssh_disconnect(ssh.session)
                                      },
                                      error = function(e){
                                          set_status(paste("Ошибка",e$message))
                                      })
                                      set_status("Завершено")
                                      enable("apply")
                                      enable("reset")
                                  })
                
                # Return something other than the promise so shiny remains responsive
                NULL
            }
            else{
                stop("Incorrect input. Sequences must be written in one-letter code.")
            }
        },
        error = function(e){
            shinyalert("Error",
                       paste("Ошибка: ",e),
                       type = "error"
            )
        })
    })
    
    output$angel.logs <- renderText({
        ### 1sec refresh
        invalidateLater(1000, session)
        validate(need(file.exists(run.logs), message = "Запустите программу."))
        get_status() %>% tail(30)
    },sep = "\n")
    
    
    output$result.proteasome = renderDataTable({
        validate(need(final.result$proteasome, message = "Запустите программу."))
        final.result$proteasome
    })
    
    output$result.TAP = renderDataTable({
        validate(need(final.result$TAP, message = "Запустите программу."))
        final.result$TAP
    })
    
    output$result.HLA = renderDataTable({
        validate(need(final.result$HLA, message = "Запустите программу."))
        final.result$HLA
    })
    
    observeEvent(input$reset,{
        session$reload()
        shinyalert("Debug",
                   as.character(file.exists(input.dst)),
                   type = "info"
        )
    })
    
    stopOnRemote = function(){
        
    } 
    
    session$onSessionEnded(function() {
        #file.remove(input.dst)
        #file.remove(run.logs)
        stopOnRemote()
        print(paste(uuid, 'the session has ended'))
    })
})
