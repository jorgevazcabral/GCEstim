process_dataset <- function(file_path, file_extension) {
  tryCatch({

    ext <- tools::file_ext(file_extension)

    data <- switch(ext,
                   csv = data.table::fread(file_path, data.table = FALSE),
                   xls = read_excel(file_path, sheet =  1),
                   xlsx = read_excel(file_path, sheet =  1),
                   validate("Invalid file; Please upload a .xls, .xlsx, .csv or .txt file")
    )

    if(ncol(data) < 3)
      stop("Data set must have at least 3 columns")

    return(data)
  }, error = function(e) {
    stop(paste("Error processing dataset:", e$message))
  })
}

server <- function(input, output, session) {

  # output$ui.data ####
  output$ui.data <- renderUI({
    if (isTRUE(input$datafile)) {
      fileInput(
        inputId = "data",
        label = "Upload dataset (XLS, XLSX, CSV, TXT)",
        accept = c("xls",".xlsx",".csv", ".txt"),
        multiple = FALSE)}
    else {
      pickerInput(
        inputId = "data",
        label = "Environment dataset",
        choices = ls(envir = globalenv()),
        options = list(style = "btn-primary"),
        multiple = FALSE,
        width = "100%")
    }
  })

  # output$ui.support.signal.points ####
  output$ui.support.signal.points <- renderUI({
    req(input$support.signal.points.n, input$support.noise.points.n, input$weight)

    if (!isTRUE(use.cv.lmgce())) {
      textInput(
        inputId = "support.signal.points",
        label = "Prior weights for the signal (comma delimited with sum equal to 1)",
        value = "0.2,0.2,0.2,0.2,0.2",
        width = "100%"
      )
    }
  })

  # output$ui.support.noise.points ####
  output$ui.support.noise.points <- renderUI({
    req(input$support.signal.points.n, input$support.noise.points.n, input$weight)

    if (!isTRUE(use.cv.lmgce())) {
      textInput(
        inputId = "support.noise.points",
        label = "Prior weights for the noise (comma delimited with sum equal to 1)",
        value = "1/3,1/3,1/3",
        width = "100%"
      )
    }
  })

  # output$ui.support.signal.pre ####
  output$ui.support.signal.pre <- renderUI({
    if (input$support.method == "ridge") {
      pickerInput(
        inputId = "support.signal.pre",
        label = "",
        choices = list("Two or more supports" = 2),
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
      )
    } else {
      pickerInput(
        inputId = "support.signal.pre",
        label = "",
        choices = list(
          "Two or more supports" = 2,
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
      )
    }
  })

  # output$ui.errormeasure.which ####
  output$ui.errormeasure.which <- renderUI({
    if (input$support.signal.pre == 2) {
      if (isTRUE(input$cv)) {
        pickerInput(
          inputId = "errormeasure.which",
          label = "Value of the error measure",
          choices = c("min", "1se", "elbow"),
          selected = "1se",
          options = list(style = "btn-primary"),
          multiple = FALSE,
          width = "100%")
      } else {
        pickerInput(
          inputId = "errormeasure.which",
          label = "Value of the error measure",
          choices = c("min", "elbow"),
          selected = "min",
          options = list(style = "btn-primary"),
          multiple = FALSE,
          width = "100%")
      }
    }
  })


  # output$ui.support.signal.equalrange ####
  output$ui.support.signal.equalrange <- renderUI(
    if (input$support.signal.pre == 1) {
      switchInput(
        inputId = "support.signal.equalrange",
        label = "Equal range",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success",
        offStatus = "danger",
        size = "mini",
        labelWidth = "150px",
        width = "250px"
      )
    }
  )

  # output$ui.support.signal.equalrange.std ####
  output$ui.support.signal.equalrange.std <- renderUI(
    if (input$support.signal.pre == 1 && isTRUE(input$support.signal.equalrange)) {
      switchInput(
        inputId = "support.signal.equalrange.std",
        label = "Standardized data",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success",
        offStatus = "danger",
        size = "mini",
        labelWidth = "150px",
        width = "250px"
      )
    }
  )

  # output$ui.support.signal.limits ####
  output$ui.support.signal.limits <- renderUI(
    if (input$support.signal.pre == 1 && isTRUE(input$support.signal.equalrange)) {
      if (isTRUE(input$support.signal.equalrange.std)) {
        numericInput(
          inputId = "support.signal.limits.std",
          label = "Upper limit of the symmetric support spaces (standardized data)",
          value = 10,
          min = 0.1,
          max = 10000,
          step = 0.1,
          width = "100%"
        )
      } else {
        numericRangeInput(
          inputId = "support.signal.limits",
          label = "Range of the support spaces (original data)",
          value = c(-10, 10),
          separator = " to ",
          min = NA,
          max = NA,
          step = NA,
          width = "100%")
      }
    } else {if (input$support.signal.pre == 1 && !isTRUE(input$support.signal.equalrange))
    {
      tags <- tagList()
      for (i in seq_len(ncol(model.matrix(as.formula(form_expr()),
                                          data = res$data)))) {
        tags[[i]] <-
          numericRangeInput(
            inputId = paste0("support.signal.limits",
                             ifelse(isTRUE(input$intercept),
                                    i-1,
                                    i)),
            label = paste0("Range of the support spaces (original data) - Beta",
                           ifelse(isTRUE(input$intercept),
                                  i-1,
                                  i)),
            value = c(-10, 10),
            separator = " to ",
            min = NA,
            max = NA,
            step = NA,
            width = "100%"
          )
      }
      tags
    }

    })


  # output$ui.support.signal.vector.given ####
  output$ui.support.signal.vector.given <- renderUI(
    if (input$support.signal.pre == 2) {
      switchInput(
        inputId = "support.signal.vector.given",
        label = "Given set",
        value = FALSE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success",
        offStatus = "danger",
        size = "mini",
        labelWidth = "150px",
        width = "250px"
      )
    }
  )

  # output$ui.support.signal.vector.n ####
  output$ui.support.signal.vector.n <- renderUI(
    if (input$support.signal.pre == 2 && !isTRUE(input$support.signal.vector.given)) {
      numericInput(
        inputId = "support.signal.vector.n",
        label = "Number of support spaces",
        value = 20,
        min = 2,
        max = 10000,
        step = 1,
        width = "100%"
      )
    }
  )

  # output$ui.support.signal.vector.limits ####
  output$ui.support.signal.vector.limits <- renderUI(
    if (input$support.signal.pre == 2 && !isTRUE(input$support.signal.vector.given)) {
      numericRangeInput(
        inputId = "support.signal.vector.limits",
        label = "Flexibility parameter",
        value = c(0.3, 20),
        separator = " to ",
        min = 0.2,
        max = NA,
        step = 0.1,
        width = "100%"
      )
    })

  # output$ui.support.signal.vector.given.which ####
  output$ui.support.signal.vector.given.which <- renderUI(
    if (input$support.signal.pre == 2 && isTRUE(input$support.signal.vector.given)) {
      textInput(
        inputId = "support.signal.vector.given.which",
        label = "UL of supports (comma delimited)",
        "0.3,0.35,0.4,0.45,0.5,0.6,0.75,1,1.25,1.5,1.75,2,3,4,5,6,8,10,12,16,20"
      )
    })

  # output$ui.support.method.ridge ####
  output$ui.support.method.ridge <- renderUI({
    if (input$support.method == "ridge") {
      tagList(
        numericInput(
          inputId = "support.method.ridge.lambda.min",
          label = "Ridge lambda min",
          value = 0.001,
          min = 0,
          step = 0.001,
          width = "100%"
        ),
        numericInput(
          inputId = "support.method.ridge.lambda.max",
          label = "Ridge lambda max",
          value = 1000,
          min = 0,
          step = 1,
          width = "100%"
        ),
        numericInput(
          inputId = "support.method.ridge.lambda.n",
          label = "Number of ridge lambdas",
          value = 100,
          min = 2,
          max = 10000,
          step = 1,
          width = "100%"
        ),
        switchInput(
          inputId = "support.method.ridge.standardize",
          label = "Standardize ridge regression",
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
          inputId = "support.method.ridge.symm",
          label = "Symmetric ridge supports",
          value = TRUE,
          onLabel = "TRUE",
          offLabel = "FALSE",
          onStatus = "success",
          offStatus = "danger",
          size = "mini",
          labelWidth = "150px",
          width = "250px"
        )
      )
    }
  })

  # output$ui.support.noise.choice ####
  output$ui.support.noise.choice <- renderUI({
    if (input$support.method == "ridge") {
      switchInput(
        inputId = "support.method.ridge.maxresid",
        label = "Maximum ridge residuals",
        value = TRUE,
        onLabel = "TRUE",
        offLabel = "FALSE",
        onStatus = "success",
        offStatus = "danger",
        size = "mini",
        labelWidth = "150px",
        width = "250px"
      )
    } else {
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
      )
    }
  })

  # output$ui.support.noise.limits ####
  output$ui.support.noise.limits <- renderUI({
    if (input$support.method == "standardized" && !isTRUE(input$support.noise.3sig)){
      if (input$support.signal.pre == 2 || (input$support.signal.pre == 1 && isTRUE(input$support.signal.equalrange) && isTRUE(input$support.signal.equalrange.std))) {
        numericRangeInput(
          inputId = "support.noise.limits",
          label = "Range of the support space (standardized data)",
          value = c(-3, 3),
          separator = " to ",
          min = NA,
          max = NA,
          step = NA,
          width = "100%")
      } else {
        numericRangeInput(
          inputId = "support.noise.limits",
          label = "Range of the support space (original data)",
          value = c(-3, 3),
          separator = " to ",
          min = NA,
          max = NA,
          step = NA,
          width = "100%")
      }
    }
  })

  # output$ui.boot.B ####
  output$ui.boot.B <- renderUI(if (isTRUE(input$boot)) {
    numericInput(
      inputId = "boot.B",
      label = "Replicates",
      value = 100,
      min = 10,
      max = 100000,
      step = 1,
      width = "100%"
    )
  })

  # output$ui.boot.method ####
  output$ui.boot.method <- renderUI(if (isTRUE(input$boot)) {
    pickerInput(
      inputId = "boot.method",
      label = "Bootstrap method",
      choices = c("residuals", "cases", "wild"),
      selected = "residuals",
      options = list(style = "btn-primary"),
      multiple = FALSE,
      width = "100%"
    )
  })

  # output$ui.cv.nfolds ####
  output$ui.cv.nfolds <- renderUI(if (isTRUE(input$cv)) {
    numericInput(
      inputId = "cv.nfolds",
      label = "Folds",
      value = 5,
      min = 3,
      max = 20,
      step = 1,
      width = "100%"
    )
  })

  # output$ui.seed ####
  output$ui.seed <- renderUI(if (isTRUE(input$cv)) {
    numericInput(
      inputId = "seed",
      label = "Seed",
      value = 230676,
      min = 1,
      max = 100000,
      step = 1,
      width = "100%"
    )
  })

  res <-
    reactiveValues(data = NULL,
                   lmgce = NULL,
                   summary = NULL,
                   expression = NULL,
                   model = NULL)

  use.cv.lmgce <- reactive({
    length(unlist(strsplit(input$support.signal.points.n, ","))) > 1 ||
      length(unlist(strsplit(input$support.noise.points.n, ","))) > 1 ||
      length(unlist(strsplit(input$weight, ","))) > 1
  })

  observeEvent(input$data, {
    withProgress(message = 'Loading dataset...', {

      if (isTRUE(input$datafile)) {
        tryCatch({
          res$data <- process_dataset(input$data$datapath,
                                      input$data$name)

          showNotification("Data set loaded successfully", type = "message")
        }, error = function(e) {
          showNotification(paste("Error:", e$message), type = "error")
          res$data <- NULL
        })
      } else {
        res$data <- try({
          dat <- get(x = input$data, envir = globalenv())
          if (inherits(dat, what = "sf")) {
            dat
          } else {
            as.data.frame(dat)
          }
        }, silent = TRUE)
      }

      updatePickerInput(session, "dep_variable",
                        choices = colnames(res$data),
                        selected = colnames(res$data)[1])

      updateActionButton(
        session,
        inputId = "n",
        label = paste0("n = ", nrow(res$data))
      )

      updateActionButton(
        session,
        inputId = "k",
        label = paste0("k = ", ncol(res$data))
      )
    })
  })

  observeEvent(input$dep_variable, {

    updatePickerInput(session, "indep_variable",
                      choices = setdiff(colnames(res$data),
                                        input$dep_variable),
                      selected = setdiff(colnames(res$data),
                                         input$dep_variable))
  })

  observeEvent(input$indep_variable, {

    updatePickerInput(session, "interaction_variable",
                      choices = input$indep_variable,
                      selected = NULL)
  })

  # output$viewdata ####
  output$viewdata <- DT::renderDataTable({
    req(res$data)
    DT::datatable(res$data,
                  options = list(
                    lengthChange = FALSE,
                    scrollX = TRUE,
                    pageLength = 20,
                    #lengthMenu = c(5, 10, 15, 20),
                    searching = FALSE,
                    info = FALSE,
                    keys = TRUE,
                    select = list(style = 'os', items = 'row'),
                    dom = 'Blfrtip',
                    rowId = 0,
                    initComplete = DT::JS(
                      "function(settings, json) {",
                      "$(this.api().table().header()).css({'background-color': '#FFA84C', 'color': 'black'});",
                      "}")
                  ),
                  autoHideNavigation = FALSE,
                  filter = "none",
                  rownames = TRUE) %>%
      DT::formatStyle(2,
                      target='row',
                      backgroundColor = "#FFE68C")
  }, serve = TRUE)

  form_expr <- reactive(paste0("`",
                               as.character(input$dep_variable),
                               "`",
                               " ~ ",
                               if (input$intercept) {"1 + `"} else {"-1 + `"},
                               paste0(
                                 as.character(input$indep_variable),
                                 collapse = "` + `"),
                               "`",
                               if (!is.null(input$interaction_variable)) {
                                 paste(
                                   "+ `",
                                   paste(as.character(input$interaction_variable),
                                         collapse = "`*`",
                                         sep = ""),
                                   "`",
                                   sep = "")
                               })
  )

  # input$run_model ####
  observeEvent(input$run_model, {
    withProgress(message = ifelse(isTRUE(use.cv.lmgce()),
                                  'Obtaining cv.lmgce model...',
                                  'Obtaining lmgce model...'), {
                                    tryCatch({
                                      session$sendCustomMessage('timer', list(id = "scenarioRuntime", event = "start"))
                                      on.exit(session$sendCustomMessage('timer', list(id = "scenarioRuntime", event = "end")))
                                      res$lmgce <-
                                        {
                                          fun_lmgce <- if (isTRUE(use.cv.lmgce())) cv.lmgce else lmgce

                                          fun_lmgce(
                                            formula = as.formula(form_expr()),
                                            data = res$data,
                                            cv = input$cv,
                                            cv.nfolds =
                                              {if (is.null(input$cv.nfolds)) 5 else input$cv.nfolds},
                                            support.method = input$support.method,
                                            support.method.ridge.lambda.min =
                                              {if (input$support.method == "ridge") input$support.method.ridge.lambda.min else NULL},
                                            support.method.ridge.lambda.max =
                                              {if (input$support.method == "ridge") input$support.method.ridge.lambda.max else NULL},
                                            support.method.ridge.lambda.n =
                                              {if (input$support.method == "ridge") input$support.method.ridge.lambda.n else NULL},
                                            support.method.ridge.standardize =
                                              {if (input$support.method == "ridge") input$support.method.ridge.standardize else NULL},
                                            support.method.ridge.symm =
                                              {if (input$support.method == "ridge") input$support.method.ridge.symm else NULL},
                                            support.method.ridge.maxresid =
                                              {if (input$support.method == "ridge") input$support.method.ridge.maxresid else NULL},
                                            support.signal = {if (input$support.signal.pre == 1) {
                                              if (isTRUE(input$support.signal.equalrange)) {
                                                if (isTRUE(input$support.signal.equalrange.std))
                                                  input$support.signal.limits.std else
                                                    input$support.signal.limits} else
                                                      matrix(
                                                        {aux.support.signal.limits <- NULL
                                                        for (i in seq_len(
                                                          ncol(
                                                            model.matrix(as.formula(form_expr()),
                                                                         data = res$data)
                                                          )
                                                        )) {
                                                          aux.support.signal.limits <-
                                                            c(aux.support.signal.limits, input[[paste0("support.signal.limits",
                                                                                                       ifelse(isTRUE(input$intercept),
                                                                                                              i-1,
                                                                                                              i))]])
                                                        }
                                                        aux.support.signal.limits
                                                        },
                                                        ncol = 2,
                                                        byrow = TRUE)
                                            } else {NULL}},
                                            support.signal.vector =
                                              {if (input$support.signal.pre == 2 && isTRUE(input$support.signal.vector.given)) {
                                                sort(as.numeric(unlist(strsplit(input$support.signal.vector.given.which, ","))),
                                                     decreasing = TRUE)
                                              }
                                                else
                                                  NULL},
                                            support.signal.vector.min =
                                              {if (is.null(input$support.signal.vector.limits)) 0.5 else
                                                input$support.signal.vector.limits[1]},
                                            support.signal.vector.max =
                                              {if (is.null(input$support.signal.vector.limits)) 20 else
                                                input$support.signal.vector.limits[2]},
                                            support.signal.vector.n =
                                              {if (is.null(input$support.signal.vector.limits)) 20 else
                                                input$support.signal.vector.n},
                                            errormeasure = input$errormeasure,
                                            errormeasure.which =
                                              {if (input$support.signal.pre == 1)
                                                "min"
                                                else
                                                  input$errormeasure.which},
                                            support.signal.points =
                                              {unlistM <- unlist(strsplit(input$support.signal.points.n, ","))
                                              q1.vect <- NULL

                                              if (isTRUE(use.cv.lmgce())) {
                                                for (i in 1:length(unlistM)){
                                                  q1.vect <- c(q1.vect,
                                                               eval(parse(text = unlistM[i])))
                                                }
                                              } else {
                                                unlistq1 <- unlist(strsplit(input$support.signal.points, ","))
                                                for (i in 1:length(unlistq1)){
                                                  q1.vect <- c(q1.vect,
                                                               eval(parse(text = unlistq1[i])))
                                                }
                                              }
                                              q1.vect
                                              },
                                            twosteps.n = input$twosteps.n,
                                            support.noise =
                                              {if (input$support.method == "ridge") NULL
                                                else if (isTRUE(input$support.noise.3sig)) NULL
                                                else
                                                  input$support.noise.limits
                                              },
                                            support.noise.points =
                                              {unlistJ <- unlist(strsplit(input$support.noise.points.n, ","))
                                              q2.vect <- NULL

                                              if (isTRUE(use.cv.lmgce())) {
                                                for (i in 1:length(unlistJ)){
                                                  q2.vect <- c(q2.vect,
                                                               eval(parse(text = unlistJ[i])))
                                                }
                                              } else {
                                                unlistq2 <- unlist(strsplit(input$support.noise.points, ","))
                                                for (i in 1:length(unlistq2)){
                                                  q2.vect <- c(q2.vect,
                                                               eval(parse(text = unlistq2[i])))
                                                }
                                              }
                                              q2.vect
                                              },
                                            weight =
                                              {unlist.weight <- unlist(strsplit(input$weight, ","))
                                              weight.vect <- NULL
                                              for (i in 1:length(unlist.weight)){
                                                weight.vect <- c(weight.vect,
                                                                 eval(parse(text = unlist.weight[i])))
                                              }
                                              weight.vect
                                              },
                                            method = input$method,
                                            caseGLM = input$caseGLM,
                                            boot.B =
                                              {if (isTRUE(input$boot)) input$boot.B else 0},
                                            boot.method =
                                              {if (isTRUE(input$boot)) input$boot.method else "residuals"},
                                            seed =
                                              {if (is.null(input$seed)) 230676 else input$seed},
                                            OLS = input$OLS)
                                        }

                                      if (isTRUE(use.cv.lmgce())) {res$model <- res$lmgce$best}
                                      else {res$model <- res$lmgce}

                                      res$summary <-
                                        summary(res$model,
                                                call = FALSE,
                                                ci.level = 0.95)

                                      if (input$support.signal.pre == 1) {
                                        updateActionButton(
                                          session,
                                          inputId = "min",
                                          label = "min",
                                          icon("circle-question",
                                               class = "fa-solid",
                                               style = "color: rgb(255,255,0)"))
                                        updateActionButton(
                                          session,
                                          inputId = "1se",
                                          label = "1se",
                                          icon("circle-question",
                                               class = "fa-solid",
                                               style = "color: rgb(255,255,0)"))
                                      }
                                      else {
                                        updateActionButton(
                                          session,
                                          inputId = "min",
                                          label = paste0("min = ",
                                                         round(res$model$support.signal.min, 2)),
                                          icon = {if (res$model$support.signal.min == max(res$model$support.ok)) {
                                            icon("thumbs-down",
                                                 class = "fa-solid",
                                                 style = "color: rgb(255,0,0)")} else {
                                                   icon("thumbs-up",
                                                        class = "fa-solid",
                                                        style = "color: rgb(0,255,0)")}}
                                        )
                                        updateActionButton(
                                          session,
                                          inputId = "1se",
                                          label = "1se",
                                          icon("circle-question",
                                               class = "fa-solid",
                                               style = "color: rgb(255,255,0)"))

                                        if (isTRUE(input$cv))
                                          updateActionButton(
                                            session,
                                            inputId = "1se",
                                            label = paste0("1se = ",
                                                           round(res$model$support.signal.1se, 2)),
                                            icon = {if (res$model$support.signal.1se == min(res$model$support.ok)) {
                                              icon("thumbs-down", class = "fa-solid", style = "color: rgb(255,0,0)")} else {
                                                icon("thumbs-up", class = "fa-solid", style = "color: rgb(0,255,0)")}}
                                          )}

                                      updateActionButton(
                                        session,
                                        inputId = "NormEnt",
                                        label = paste0("NormEnt = ",
                                                       round(res$model$nep,3),
                                                       " (",
                                                       ifelse(input$cv,
                                                              round(res$model$nep.cv.mean,3),
                                                              "-"),
                                                       ")"),
                                        icon = {if (res$model$convergence == 1) {
                                          icon("thumbs-down", class = "fa-solid", style = "color: rgb(255,0,0)")} else {
                                            icon("thumbs-up", class = "fa-solid", style = "color: rgb(0,255,0)")}}
                                      )

                                      updateActionButton(
                                        session,
                                        inputId = "R2",
                                        label = paste0("R2 = ",
                                                       round(res$summary$r.squared,3)),
                                        icon = {if (res$model$convergence == 1) {
                                          icon("thumbs-down", class = "fa-solid", style = "color: rgb(255,0,0)")} else {
                                            icon("thumbs-up", class = "fa-solid", style = "color: rgb(0,255,0)")}}
                                      )

                                      updateActionButton(
                                        session,
                                        inputId = "Error",
                                        label = paste0(res$model$error,
                                                       " = ",
                                                       round(res$model$error.measure,3),
                                                       " (",
                                                       ifelse(input$cv,
                                                              round(res$model$error.measure.cv.mean,3),
                                                              "-"),
                                                       ")"),
                                        icon = {if (res$model$convergence == 1) {
                                          icon("thumbs-down", class = "fa-solid", style = "color: rgb(255,0,0)")} else {
                                            icon("thumbs-up", class = "fa-solid", style = "color: rgb(0,255,0)")}}
                                      )

                                    }, error = function(e) {
                                      showNotification(paste("Error:", e$message),
                                                       duration = NULL,
                                                       type = "error")
                                    })
                                  })
  })

  # output$plot ####
  observeEvent(c(input$run_model,
                 input$ci.method,
                 input$ci.level,
                 input$OLS.plot), {

                   output$plot1 <- renderPlot({
                     req(res$model)
                     plot(res$model,type = "ggplot2",
                          which = 1,
                          ci.level = {if (isTRUE(input$boot)) input$ci.level else 0.95},
                          ci.method = {if (isTRUE(input$boot)) input$ci.method else "z"},
                          OLS = input$OLS.plot)
                   })
                 })

  observeEvent(c(input$run_model,
                 input$NormEnt.plot,
                 input$OLS.plot), {
                   output$plot2 <- renderPlot({
                     req(res$model)
                     plot(res$model,type = "ggplot2",
                          which = 2,
                          OLS = input$OLS.plot,
                          NormEnt = input$NormEnt.plot)
                   })
                 })

  observeEvent(c(input$run_model,
                 input$OLS.plot), {
                   output$plot3 <- renderPlot({
                     req(res$model)
                     plot(res$model,type = "ggplot2",
                          which = 3,
                          OLS = input$OLS.plot)
                   })

                   output$plot6 <- renderPlot({
                     req(res$model)
                     plot(res$model,type = "ggplot2",
                          which = 6,
                          OLS = input$OLS.plot)
                   })
                 })

  output$plot4 <- renderPlot({
    req(res$model)
    plot(res$model,type = "ggplot2",
         which = 4)
  })

  # output$summary ####
  output$summary <- renderPrint({
    req(res$summary)
    if(!is.null(res$summary))
      res$summary
  })

  res$expression <- reactive(
    paste0(
      if (isTRUE(use.cv.lmgce()))
        "GCEstim::cv.lmgce(formula = "
      else
        "GCEstim::lmgce(formula = ",
      form_expr(),
      ",\n data = ",
      {if (!isTRUE(input$datafile)) input$data
        else
          input$data$name},
      ", cv = ",
      input$cv,
      {if (!is.null(input$cv.nfolds) && isTRUE(input$cv))
        paste0(
          ", cv.nfolds = ",
          input$cv.nfolds)},
      ", support.method = \"",
      input$support.method,
      "\"",
      ", errormeasure = \"",
      input$errormeasure,
      {if (input$support.signal.pre == 2) {
        if (is.null(input$errormeasure.which))
          paste0("\", errormeasure.which = \"1se")
        else
          paste0(
            "\", errormeasure.which = \"",
            input$errormeasure.which)}},
      "\"",
      {if (input$support.signal.pre == 1) {
        if (isTRUE(input$support.signal.equalrange)) {
          if (isTRUE(input$support.signal.equalrange)) {
            paste0(", support.signal = ",
                   input$support.signal.limits.std)
          } else {
            paste0(
              ", support.signal = c(",
              paste0(input$support.signal.limits[1],
                     ", ",
                     input$support.signal.limits[2]),
              ")")}} else {
                paste0(
                  ", support.signal = matrix(c(",
                  {aux.support.signal.limits <- NULL
                  for (i in seq_len(
                    ncol(
                      model.matrix(as.formula(form_expr()),
                                   data = res$data)
                    )
                  )) {
                    aux.support.signal.limits <-
                      c(aux.support.signal.limits, input[[paste0("support.signal.limits",
                                                                 ifelse(isTRUE(input$intercept),
                                                                        i-1,
                                                                        i))]])
                  }
                  paste(aux.support.signal.limits, collapse = ", ")
                  },
                  "), ncol = 2, byrow = TRUE)")
              }
      } else {
        paste0(", support.signal = NULL")
      }
      },
      {if (!is.null(input$support.signal.vector.limits) && input$support.signal.pre == 2) {
        if (!isTRUE(input$support.signal.vector.given)) {
          paste0(
            ", support.signal.vector.min = ",
            input$support.signal.vector.limits[1],
            ", support.signal.vector.max = ",
            input$support.signal.vector.limits[2],
            ", support.signal.vector.n = ",
            input$support.signal.vector.n)}
      }},
      {if (input$support.signal.pre == 2 && isTRUE(input$support.signal.vector.given))
      {paste0(", support.signal.vector = c(",
              input$support.signal.vector.given.which,
              ")"
      )
      }
      },
      ", support.signal.points = c(",
      {if (isTRUE(use.cv.lmgce()))
        input$support.signal.points.n
        else
          input$support.signal.points},
      ")",
      ", support.noise = ",
      {if (input$support.method == "ridge") paste0("NULL")
        else if (isTRUE(input$support.noise.3sig)) paste0("NULL")
        else
          paste0("c(",
                 input$support.noise.limits[1],
                 ", ",
                 input$support.noise.limits[2],
                 ")")
      },
      ", support.noise.points = c(",
      {if (isTRUE(use.cv.lmgce()))
        input$support.noise.points.n
        else
          input$support.noise.points},
      "), weight = c(",
      input$weight,
      ")",
      ", twosteps.n = ",
      input$twosteps.n,
      ", caseGLM = \"",
      input$caseGLM,
      "\", method = \"",
      input$method,
      "\"",
      {if (input$support.method == "ridge")
        paste0(
          ", support.method.ridge.lambda.min = ",
          input$support.method.ridge.lambda.min,
          ", support.method.ridge.lambda.max = ",
          input$support.method.ridge.lambda.max,
          ", support.method.ridge.lambda.n = ",
          input$support.method.ridge.lambda.n,
          ", support.method.ridge.standardize = ",
          input$support.method.ridge.standardize,
          ", support.method.ridge.symm = ",
          input$support.method.ridge.symm,
          ", support.method.ridge.maxresid = ",
          input$support.method.ridge.maxresid
        )
      },
      {if (isTRUE(input$boot)) paste0(", boot.B = ",
                                      input$boot.B,
                                      ", boot.method = \"",
                                      input$boot.method,
                                      "\"")},
      {if (!is.null(input$seed) && isTRUE(input$cv))
        paste0(
          ", seed = ",
          input$seed)},
      ", OLS = ",
      input$OLS,
      ")"))

  # output$code ####
  output$code <- renderUI({
    req(res$expression())
    HTML(highlight(res$expression(),
                   code = TRUE,
                   classes = classes_pandoc()))
  })

  observeEvent(input$done, {
    rstudioapi::insertText(res$expression())
    stopApp()
  })
}
