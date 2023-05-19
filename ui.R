#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(shinyTree)

ui <- tagList(
    useShinyjs(),
    navbarPage("ANGEL",
               #rmarkdown::render("intro_rus.Rmd", rmarkdown::html_document(toc = TRUE,toc_float = T))
               tabPanel("Введение", fluidPage(htmlOutput("intro"))),
               tabPanel("Сервис",
                        fluidPage(
                            sidebarLayout(
                                sidebarPanel(
                                    wellPanel(
                                        textAreaInput("input.text",label = "Введите последовательности (1 на строку)"),
                                    ),
                                    wellPanel(
                                        sliderInput("proteasome.threshold", label="Порог для протеасомы", min = 0, max = 1, step = 0.01, value = 0.5),#TODO - change value to optimal
                                        sliderInput("TAP.threshold", label="Порог для TAP", min = 0, max = 1, step = 0.01, value = 0.5),#TODO - change value to optimal
                                        sliderInput("HLA.threshold", label="Порог для HLA", min = 0, max = 1, step = 0.01, value = 0.5),#TODO - change value to optimal
                                    ),
                                    checkboxInput("include.TAP",label = "Следует ли учитывать прогноз для TAP?", value = T),
                                    wellPanel(
                                        fluidRow(
                                            column(4, align = "center", actionButton("apply",label = "Запустить")),
                                            column(4, align = "center", actionButton("reset",label = "Очистить")),
                                            column(4, align = "center", disabled(downloadButton("download", label = "Скачать")))
                                        )
                                    )
                                ),
                                mainPanel(
                                    shinyTree("result")
                                )
                            )
                        ))
    )
)
