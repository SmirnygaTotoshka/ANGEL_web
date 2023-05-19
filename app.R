library(shiny)
Sys.setlocale(category = "LC_ALL", locale = "ru_RU.UTF-8")
options(shiny.maxRequestSize=1*1024^2)
options(shiny.error = browser)
source("ui.R", encoding = "UTF-8")
source("server.R", encoding = "UTF-8")

shinyApp(ui = ui, server = server)