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
    uuid = UUIDgenerate(output = "string")
    input.dst = file.path(server_config["userdata_path","V2"],paste0(uuid,".csv"))
    run.logs = file.path(server_config["userdata_path","V2"],"logs",paste0(uuid,".log"))
    track_usage(storage_mode = store_json(path = file.path(server_config["userdata_path","V2"],"shiny")))
    
    
    final.result = reactiveValues(
        proteasome = NULL,
        HLA = NULL
    )


    write.input.to.file = function(sequences){
        fileConn<-file(input.dst)
        writeLines(c("proteins",sequences), fileConn)
        close(fileConn)
    }
    
    check.input = function(sequences){
        all(grepl("^[ACDEFGHIKLMNPQRSTVWY]+$",sequences))#TODO check min and max len
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
            showNotification(success_status)
        }
        else{
            showNotification("Error")
            set_status("Error")
            set_status(rawToChar(out$stdout))
            set_status(rawToChar(out$stderr))
            stop("Error. See log.")
        }
    }
    #TODO - check success PASS launch 
    pass.launch = function(ssh, config, success_status){
        launch.cmd.with.waiting(ssh.session,command = paste(server_config["pass_program","V2"], proteasome.config.path),
                                success_status = "Протеасома. Запуск.")
    }
    
    observeEvent(input$apply,{
        tryCatch({
            proteasome_threshold = isolate(input$proteasome.threshold)
            TAP_threshold = isolate(input$TAP.threshold)
            HLA_threshold = isolate(input$HLA.threshold)
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
                    shinyalert(title = paste("Запуск",uuid),text = "Задание начато к исполнению. Это займет некоторое время. Не перезагружайте страницу!", type = "info")
                    workdir = paste0(server_config["pass_root","V2"],"\\", uuid)
                    ssh_exec_wait(ssh.session, paste("mkdir",workdir))
                    #scp_upload(ssh.session, input.dst, to = paste0(uuid,".csv"))#WARNING Cannot directly load to needed folder
                    launch.cmd.with.waiting(ssh.session, command = paste0("python ",server_config["python_scripts","V2"],"\\transferFiles.py get_input ",
                                                                          paste(uuid, server_config["ftp_server","V2"], server_config["ftp_server_user","V2"], server_config["ftp_server_passwd","V2"])),
                                            success_status = "Входные последовательности готовы")
                    #ssh_exec_wait(ssh.session, paste("move",paste0("'",uuid,".csv'"), workdir))
                    ssh_exec_wait(ssh.session, paste("cd",workdir))
                    # launch.cmd.with.waiting(ssh.session,command = paste0("python ",server_config["python_scripts","V2"],"\\renameInput.py ",workdir),
                    #                         success_status = "Входные последовательности гот")#Убрать эти кавычки странные
                    remote.input.path = paste0(workdir,"\\","ready_Proteasome_", uuid, ".csv")
                    path_to_config = paste0(workdir,"\\", uuid, ".json")
                    launch.cmd.with.waiting(ssh.session, command = paste0("python ",server_config["python_scripts","V2"],"\\prepareProteasomeInput.py ",uuid),
                                            success_status = "Входной файл для моделирования готов.")                    
                    launch.cmd.with.waiting(ssh.session,command = paste0("python ",server_config["python_scripts","V2"],"\\generateConfigForConverter.py -i ",
                                                                         remote.input.path," -o ", workdir," -c peptide -u -t 5 -f ",uuid," ",path_to_config),
                                            success_status = "Есть конфигурация для конвертера.")
                    launch.cmd.with.waiting(ssh.session,command = paste0("python ",server_config["python_scripts","V2"],"\\SeqToSDF.py ",path_to_config),
                                            success_status = "Сконвертировано в SDF.")
                    #invivo windows_radius =3 MNA level = 9
                    proteasome.config.path = paste0(server_config["pass_root","V2"], "\\", uuid,"_prot_pred.txt")
                    launch.cmd.with.waiting(ssh.session,command = paste0("python ",server_config["python_scripts","V2"],"\\generateConfigForPrediction.py -m models/Proteasome.MSAR -s ",
                                                                         uuid,"/", uuid, ".sdf"," -a id -p ",uuid," ", proteasome.config.path),
                                            success_status = "Протеасома. Готова к запуску.")
                    launch.cmd.with.waiting(ssh.session,command = paste(server_config["pass_program","V2"], proteasome.config.path),
                                            success_status = "Протеасома. Запуск завершен.")
                    launch.cmd.with.waiting(ssh.session, command = paste0("python ",server_config["python_scripts","V2"],"\\processProteasomeOutput.py ",
                                                                          paste(uuid, proteasome_threshold)),
                                            success_status = "Результат протеасомы обработан.")
                    launch.cmd.with.waiting(ssh.session, command = paste0("python ",server_config["python_scripts","V2"],"\\transferFiles.py send_proteasome ",
                                                                          paste(uuid, server_config["ftp_server","V2"], server_config["ftp_server_user","V2"], server_config["ftp_server_passwd","V2"])),
                                            success_status = "Результат протеасомы получен.")
                    #host = str_split(pass_server,"\\@")[[1]][2]
                    #user = str_split(pass_server,"\\@")[[1]][1]
                    #f = RCurl::scp(host = host, path = paste0(workdir,"\\", uuid, "_Proteasome_result.csv"), keypasswd = pass_passwd,user = user,binary =  F)
                    #scp_download(ssh.session, files = paste0(workdir,"\\", uuid, "_Proteasome_result.csv"), to = server_config["userdata_path","V2"])
                    #scp_download(ssh.session, files = paste0("C:\\Users\\SmirnygaTotoshka\\Desktop\\diplom\\angel_userdata\\PASS2CSV.log"), to = server_config["userdata_path","V2"])
                    final.result$proteasome = read.csv2(file.path(server_config["userdata_path","V2"],paste0(uuid,"_Proteasome.csv")))
                    final.result$HLA = data.frame(matrix(1:10,nrow = 5))
                })
                
                # Catch inturrupt (or any other error) and notify user
                result <- catch(result,
                                function(e){
                                    final.result$proteasome = NULL
                                    final.result$HLA = NULL
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
            enable("apply")
            enable("reset")
        })
    })
    
    output$angel.logs <- renderText({
        ### 1sec refresh
        invalidateLater(100, session)
        validate(need(file.exists(run.logs), message = "Запустите программу."))
        get_status() %>% tail(30)
    },sep = "\n")
    
    
    output$result.proteasome = renderDataTable({
        validate(need(final.result$proteasome, message = "Запустите программу."))
        final.result$proteasome
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
