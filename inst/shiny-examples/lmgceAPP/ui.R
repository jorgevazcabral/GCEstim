# library(shiny)
# library(miniUI)
# library(shinyWidgets)
# library(shinydashboard)
# library(shinydashboardPlus)
# library(readxl)
# library(DT)
# library(ggplot2)
# library(plotly)
# library(downlit)

ui <- miniPage(
  gadgetTitleBar(
    "Generalized Cross Entropy linear regression",
    right = miniTitleBarButton("done", "Export", primary = TRUE)),
  miniButtonBlock(
    actionButton(
      inputId = "run_model",
      label = "Run model",
      width = "100%",
      class = "btn-success btn-sm"
    ),
    actionButton(inputId = "min",
                 label = "min = ",
                 icon("circle-question", class = "fa-solid", style = "color: rgb(255,255,0)"),
                 disabled = TRUE,
                 style = "color: black; background-color: #7799FF; border-color: #2e6da4"),
    actionButton(inputId = "1se",
                 label ="1se = ",
                 icon("circle-question", class = "fa-solid", style = "color: rgb(255,255,0)"),
                 disabled = TRUE,
                 style = "color: black; background-color: #7799FF; border-color: #2e6da4"),
    p(id = "scenarioRuntime",
      tags$label(class = "minutes"),
      tags$label(class = "seconds")),
    tags$script(HTML(
      '
        $(function(){
            var timer;

            Shiny.addCustomMessageHandler("timer", function(data){
                if(data.event === "end") return clearInterval(timer);

                var minutesLabel = document.querySelector(`#${data.id} .minutes`);
                var secondsLabel = document.querySelector(`#${data.id} .seconds`);
                var totalSeconds = 0;

                function pad(val) {
                  var valString = val + "";
                  if (valString.length < 2) {
                    return "0" + valString;
                  } else {
                    return valString;
                  }
                }
                function setTime() {
                  ++totalSeconds;
                  secondsLabel.innerHTML = pad(totalSeconds % 60);
                  minutesLabel.innerHTML = `${pad(parseInt(totalSeconds / 60))} : `;
                }

                timer = setInterval(setTime, 1000);
            });
        });
        '
    ))
  ),
  miniButtonBlock(
    actionButton(
      inputId = "NormEnt",
      label = "NormEnt = ",
      icon("circle-question", class = "fa-solid", style = "color: rgb(255,255,0)"),
      disabled = TRUE,
      style = "color: black; background-color: #7799FF; border-color: #2e6da4"),
    actionButton(
      inputId = "R2",
      label = "R2 =",
      icon("circle-question", class = "fa-solid", style = "color: rgb(255,255,0)"),
      disabled = TRUE,
      style = "color: black; background-color: #7799FF; border-color: #2e6da4"),
    actionButton(
      inputId = "Error",
      label = "Error =",
      icon("circle-question", class = "fa-solid", style = "color: rgb(255,255,0)"),
      disabled = TRUE,
      style = "color: black; background-color: #7799FF; border-color: #2e6da4")
  ),
  miniTabstripPanel(
    miniTabPanel(
      "Data",
      icon = icon("database"),
      miniContentPanel(
        switchInput(
          inputId = "datafile",
          label = "file",
          value = FALSE,
          onLabel = "TRUE",
          offLabel = "FALSE",
          onStatus = "success",
          offStatus = "danger",
          size = "mini",
          labelWidth = "150px",
          width = "250px"
        ),
        uiOutput("ui.data"),
        DTOutput("viewdata"),
        miniButtonBlock(
          actionButton(
            inputId = "n",
            label = "n = ",
            disabled = TRUE,
            style = "color: black; background-color: #7799FF; border-color: #2e6da4"),
          actionButton(
            inputId = "k",
            label = "k = ",
            disabled = TRUE,
            style = "color: black; background-color: #7799FF; border-color: #2e6da4")
        ))
    ),
    miniTabPanel(
      "lmgce()",
      icon = icon("sliders"),
      miniTabstripPanel(
        miniTabPanel(
          "Formula",
          icon = icon("square-root-variable"),
          miniContentPanel(
            pickerInput(
              inputId = "dep_variable",
              label = "Dependent Variable",
              choices = NULL,
              selected = NULL,
              options = list(style = "btn-primary"),
              multiple = FALSE,
              width = "100%"
            ),
            switchInput(
              inputId = "intercept",
              label = "Intercept",
              value = TRUE,
              onLabel = "TRUE",
              offLabel = "FALSE",
              onStatus = "success",
              offStatus = "danger",
              size = "mini",
              labelWidth = "150px",
              width = "250px"
            ),
            pickerInput(
              inputId = "indep_variable",
              label = "Independent Variables",
              choices = NULL,
              selected = NULL,
              options = list(
                `actions-box` = TRUE,
                `deselect-all-text` = "None",
                `select-all-text` = "All",
                `none-selected-text` = "None",
                style = "btn-primary"
              ),
              multiple = TRUE,
              width = "100%"
            ),
            pickerInput(
              inputId = "interaction_variable",
              label = "Interaction Variables",
              choices = NULL,
              selected = NULL,
              options = list(
                `actions-box` = TRUE,
                `deselect-all-text` = "None",
                `select-all-text` = "All",
                `none-selected-text` = "None",
                style = "btn-primary"
              ),
              multiple = TRUE,
              width = "100%"
            ),
            pickerInput(
              inputId = "caseGLM",
              label = "general linear model case",
              choices = list("data" = "D",
                             "moment" = "M",
                             "normed-moment" = "NM"),
              selected = "data",
              options = list(style = "btn-primary"),
              multiple = FALSE,
              width = "100%"
            )
          )
        ),
        miniTabPanel(
          "Signal support",
          icon = icon("signal", class = "fa-solid"),
          miniContentPanel(
            textInput(
              inputId = "support.signal.points",
              label = "Prior weights for the signal (comma delimited with sum equal to 1)",
              "0.2,0.2,0.2,0.2,0.2",
              width = "100%"),
            pickerInput(
              inputId = "support.signal.pre",
              label = "",
              choices = list(
                "Two or more supports (standardized data)" = 2,
                "One support" = 1
              ),
              selected = 2,
              options = list(
                `actions-box` = TRUE,
                `deselect-all-text` = "None",
                `select-all-text` = "All",
                `none-selected-text` = "None",
                style = "btn-primary"
              ),
              multiple = FALSE,
              width = "100%"
            ),
            uiOutput("ui.errormeasure.which"),
            uiOutput("ui.support.signal.equalrange"),
            uiOutput("ui.support.signal.equalrange.std"),
            uiOutput("ui.support.signal.limits"),
            uiOutput("ui.support.signal.vector.given"),
            uiOutput("ui.support.signal.vector.given.which"),
            uiOutput("ui.support.signal.vector.n"),
            uiOutput("ui.support.signal.vector.limits")
          )
        ),
        miniTabPanel(
          "Noise support",
          icon = icon("wave-square", class = "fa-solid"),
          miniContentPanel(
            sliderInput(
              inputId = "weight",
              label = "Noise weight",
              min = 0,
              max = 1,
              value = 0.5,
              step =0.01,
              ticks = FALSE,
              width = "100%"),
            textInput(
              inputId = "support.noise.points",
              label = "Prior weights for the noise (comma delimited with sum equal to 1)",
              "1/3,1/3,1/3",
              width = "100%"),
            switchInput(
              inputId = "support.noise.3sig",
              label = "3 sigma rule",
              value = TRUE,
              onLabel = "TRUE",
              offLabel = "FALSE",
              onStatus = "success",
              offStatus = "danger",
              size = "mini",
              labelWidth = "150px",
              width = "250px"
            ),
            uiOutput("ui.support.noise.limits")
          )
        ),
        miniTabPanel(
          "Other",
          icon = icon("ellipsis"),
          miniContentPanel(
            pickerInput(
              inputId = "errormeasure",
              label = "Error measure",
              choices = c("RMSE", "MSE", "MAD","MAE", "MAPE", "sMAPE", "MASE"),
              selected = "RMSE",
              options = list(style = "btn-primary"),
              multiple = FALSE,
              width = "100%"
            ),
            numericInput(
              inputId = "twosteps.n",
              label = "Number of post GCE reestimations",
              value = 1,
              min = 0,
              max = 10000,
              step = 1,
              width = "100%"
            ),
            pickerInput(
              inputId = "method",
              label = "Optimization method",
              choices = c("dual.BFGS", "dual.CG", "dual.L-BFGS-B",
                          "dual.Rcgmin", "dual.bobyqa", "dual.newuoa",
                          "dual.nlminb", "dual.nlm",
                          "dual.lbfgs",
                          "dual.lbfgsb3c",
                          "dual",
                          "primal.solnl", "primal.solnp"),
              selected = "dual.BFGS",
              options = list(style = "btn-primary"),
              multiple = FALSE,
              width = "100%"
            ),
            switchInput(
              inputId = "OLS",
              label = "OLS",
              value = TRUE,
              onLabel = "TRUE",
              offLabel = "FALSE",
              onStatus = "success",
              offStatus = "danger",
              size = "mini",
              labelWidth = "150px",
              width = "250px"
            ),
            switchInput(
              inputId = "boot",
              label = "Bootstrapp",
              value = FALSE,
              onLabel = "TRUE",
              offLabel = "FALSE",
              onStatus = "success",
              offStatus = "danger",
              size = "mini",
              labelWidth = "150px",
              width = "250px"
            ),
            uiOutput("ui.boot.B"),
            uiOutput("ui.boot.method"),
            switchInput(
              inputId = "cv",
              label = "Cross-validation",
              value = TRUE,
              onLabel = "TRUE",
              offLabel = "FALSE",
              onStatus = "success",
              offStatus = "danger",
              size = "mini",
              labelWidth = "150px",
              width = "250px"
            ),
            #fillRow(
            uiOutput("ui.cv.nfolds"),
            uiOutput("ui.seed")
            #)
          ))
      )),
    miniTabPanel(
      "summary()",
      icon = icon("table"),
      miniContentPanel(
        verbatimTextOutput("summary")
      )
    ),
    miniTabPanel(
      "plot()",
      icon = icon("area-chart"),
      dropdownButton(
        inputId = "settings.plot",
        label = "Settings",
        icon = icon("sliders"),
        status = "primary",
        circle = FALSE,
        width = "100%",
        radioGroupButtons(
          inputId = "ci.method",
          label = "CI method",
          choices = c("z", "percentile", "basic"),
          selected = "z",
          status = "primary",
          size = "normal",
          justified = TRUE,
          width = "100%"),
        # prettyRadioButtons(
        #   inputId = "ci.method",
        #   label = "CI method",
        #   choices = c("z", "percentile", "basic"),
        #   selected = "z",
        #   status = "primary",
        #   icon = icon("check"),
        #   animation = "pulse",
        #   width = "100%",
        #   inline = TRUE,
        #   ),
        numericInput(
          inputId = "ci.level",
          label = "Confidence level",
          value = 0.95,
          min = 0.1,
          max = 0.9999,
          step = 0.0001,
          width = "100%"
        ),
        switchInput(
          inputId = "OLS.plot",
          label = "plot OLS",
          value = TRUE,
          onLabel = "TRUE",
          offLabel = "FALSE",
          onStatus = "success",
          offStatus = "danger",
          size = "mini",
          labelWidth = "150px",
          width = "250px"
        ),
        switchInput(
          inputId = "NormEnt.plot",
          label = "plot NormEnt",
          value = TRUE,
          onLabel = "TRUE",
          offLabel = "FALSE",
          onStatus = "success",
          offStatus = "danger",
          size = "mini",
          labelWidth = "150px",
          width = "250px"
        )
      ),
      miniTabstripPanel(
        miniTabPanel(
          "Estimates",
          miniContentPanel(
            plotOutput("plot1", height = "100%")
          )
        ),
        miniTabPanel(
          "Error vs support",
          miniContentPanel(
            plotOutput("plot2", height = "100%")
          )
        ),
        miniTabPanel(
          "Estimates vs support",
          miniContentPanel(
            plotOutput("plot3", height = "100%")
          )
        ),
        miniTabPanel(
          "NormEnt vs support",
          miniContentPanel(
            plotOutput("plot4", height = "100%")
          )
        ),
        miniTabPanel(
          "GCE Reestimation",
          miniContentPanel(
            plotOutput("plot6", height = "100%")
          )
        )
      )
    ),
    miniTabPanel(
      "Code",
      icon = icon("code"),
      miniContentPanel(
        uiOutput("code")
      )
    ),
    miniTabPanel(
      "About",
      icon = icon("circle-info"),
      miniContentPanel(
        h3(HTML("<b>Version 0.2 developed by:</b>")),
        h4("Jorge Cabral"),
        tags$h6(tags$a(href = "mailto:jorgecabral@ua.pt","jorgecabral@ua.pt")),
        h4("Pedro Macedo"),
        h4("Vera Afreixo"),
        hr(),
        tags$head(
          tags$style(HTML("hr {border-top: 1px solid #000000;}"))
        ),
        h6("This work is supported by CIDMA under the FCT
              (Portuguese Foundation for Science and Technology)
                 Multi-Annual Financing Program for R&D Units.")
      )
    )
  )
)
