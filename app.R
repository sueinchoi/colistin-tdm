## Colistin Bayesian Dosing App — Modern UI Redesign
## Based on the **Run330** popPK model. Lineage: Run240 -> Run241 corrected the inverted
## Cockcroft-Gault sex factor (SEX 1 = Male); Run241 -> Run242 restored ID 3's BSA to the CRF
## original (1.24; the dataset had 1.49 carried over from ID 2); Run242 -> Run249 replaced BSA
## on CL_CMS with fixed allometric weight (BSA never survives the covariate search once the sex
## factor is corrected); Run249 -> Run330 added a peripheral compartment to colistin B.
## Supersedes Run249/Run242/Run241/Run240/Run210.
##
## Two things changed here that alter predictions and must not be reverted piecemeal:
##   - BSA is no longer a covariate. Weight enters as (WT/70)^0.75 on clearances and (WT/70)
##     on volumes, so **WT must be passed to param()** — omitting it silently simulates a 70 kg
##     patient regardless of the sidebar weight.
##   - Colistin B is 2-compartment (Q_MB = 0.0894 L/h), so the terminal colistin B profile is
##     no longer mono-exponential.

suppressMessages({
  library(shiny); library(mrgsolve); library(dplyr)
  library(ggplot2); library(DT); library(tidyr)
  library(bslib); library(shinycssloaders); library(shinyWidgets)
})

# ---- Model loading ----
## Model file sits next to the app. An absolute path breaks both a fresh clone and
## any shinyapps.io deployment, because that path does not exist there.
MOD_CODE_FILE <- "colistin_run330_mrgsolve.rds"
mod_code <- readRDS(MOD_CODE_FILE)
mod <- mcode("colistin_app", mod_code, quiet = TRUE)
pop_par <- as.list(param(mod))
# Run330 OMEGA variances: CL 0.0520812, V1 0.747603, CLMA 0.466954, CLMB 1.83435
# (EDA/115_run330_mrgsolve.R prints these — regenerate and paste after any new final run.)
OMEGA_SD <- c(sqrt(0.0520812), sqrt(0.747603), sqrt(0.466954), sqrt(1.83435))

## ---- 치료 창 (2026-08-06 결합보정 기준으로 교체 · 2026-08-07 f_u 를 Run330 값으로) ----
## 이전 판은 총약물 2-8 mg/L 를 썼다. 그 8 mg/L 는 인용된 문헌(Sorli 2013, Forrest 2017)
## 어디에도 없는 값이고, 이 코호트의 측정 비결합분율은 consensus 가정 0.5 의 절반 이하다.
##
## f_u,eff 주의 — 2026-08-07 에 바로잡았다.
##   f_u,eff 는 F_MA/CL_MA : F_MB/CL_MB 로 가중되므로 **모형이 바뀌면 같이 바뀐다.**
##   Run242 에서 0.2232(-> 0.22), Run330 에서 **0.2251** 이다. 앱과 원고가 오래
##   Run242 값을 쓰고 있었다. 그리고 0.2251 은 소수 둘째자리 반올림 경계에 정확히
##   걸려 있어 접은 뒤 유도하면 목표가 4.73/4.53 으로 4% 튄다 — **접지 않고 유도한다.**
##
##   효능  fAUC24/MIC >= 25, MIC 1 mg/L, f_u 0.2251  -> 총약물 25/24/0.2251 = 4.63 mg/L
##   안전  consensus 상한 4 mg/L 를 그 문헌 코호트의 결합으로 환산
##         = 4 x f_u,ref / 0.2251. f_u,ref 는 보고된 적이 없어 민감도로 다룬다.
##         **주 가정 0.31 -> 5.51 mg/L**, 범위 0.22-0.50 -> 3.91-8.88 mg/L.
## 창이 4.63-5.51 로 좁고, 가장 보수적인 f_u,ref = 0.22 에서는 상한이 하한보다 낮아진다.
## 그 사실을 UI 에 드러낸다 — 숨기면 도구가 없는 정밀도를 주장하게 된다.
FU_EFF       <- 0.2251 # median per-patient f_u,eff, Run330 (Supplementary Table S2c)
TARGET_LOW   <- 25/24/FU_EFF      # 효능 하한 4.63 (총약물)
TARGET_HIGH  <- 4*0.31/FU_EFF     # 안전 상한 5.51 (총약물, f_u,ref = 0.31)
TARGET_HI_LO <- 4*0.22/FU_EFF     # 안전 상한 3.91, 가장 보수적인 f_u,ref = 0.22
TARGET_HI_HI <- 4*0.50/FU_EFF     # 안전 상한 8.88, f_u,ref = 0.50 (consensus 가정)

## Run330 SIGMA, as SDs on the ng/mL scale (the .ext gives variances).
## Used as the MAP likelihood's residual model — CMS / colistin A / colistin B.
RES_ADD <- c(sqrt(2.22606e+06), sqrt(96230.9), sqrt(353.609))   # 1492, 310.2, 18.8 ng/mL
RES_PRP <- c(sqrt(0.0986821),   sqrt(0.0270764), sqrt(0.0717999))  # 0.314, 0.165, 0.268

## Covariate exponents, read from the loaded model so they can never drift from the .rds.
EXP_CL   <- pop_par$EXP_CLCR_CL   # 0.367
EXP_MA   <- pop_par$EXP_CLCR_MA   # 0.551
EXP_MB   <- pop_par$EXP_CLCR_MB   # 1.149
EXP_WT_C <- pop_par$EXP_WT_CL     # 0.75 FIXED
EXP_WT_V <- pop_par$EXP_WT_V      # 1.0  FIXED

# ==================== UI ====================
ui <- page_sidebar(
  title = tags$div(
    style = "display:flex; align-items:center; gap:14px; padding:4px 0;",
    tags$i(class = "fas fa-pills",
           style = "color:#7FDBFF; font-size:32px; text-shadow: 0 1px 3px rgba(0,0,0,0.3);"),
    tags$span("Colistin Precision Dosing",
              style = "font-weight:700; color:#7FDBFF; font-size:24px; letter-spacing:0.5px; text-shadow: 0 1px 3px rgba(0,0,0,0.35);"),
    tags$span("v1.0",
              style = "font-size:13px; color:#b8d8e8; margin-left:6px; font-weight:500;")
  ),
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2E86AB",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter"),
    font_scale = 0.95
  ),

  tags$head(
    tags$link(rel = "stylesheet",
              href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.5.1/css/all.min.css"),
    tags$style(HTML("
      /* Navbar title area */
      .navbar-brand, .navbar .navbar-brand { padding: 8px 0 !important; }
      .navbar { min-height: 64px; }
      /* Global spacing */
      .card { box-shadow: 0 2px 8px rgba(0,0,0,0.06); border: none !important; border-radius: 10px; }
      .card-header { background: #f8f9fa !important; border-bottom: 2px solid #2E86AB !important; font-weight: 600; }
      /* Section headers in sidebar */
      .section-header { font-weight: 600; color: #2E86AB; margin-top: 14px; margin-bottom: 8px; padding-bottom: 4px; border-bottom: 2px solid #e9ecef; display: flex; align-items: center; gap: 8px; }
      .section-header i { color: #2E86AB; }
      /* Metric boxes */
      .metric-box { background: linear-gradient(135deg, #f8f9fa 0%, #ffffff 100%); border-left: 4px solid #2E86AB; padding: 12px 16px; border-radius: 6px; margin-bottom: 10px; }
      .metric-label { font-size: 12px; color: #6c757d; text-transform: uppercase; letter-spacing: 0.5px; }
      .metric-value { font-size: 28px; font-weight: 700; color: #2E86AB; line-height: 1.1; }
      .metric-sub { font-size: 13px; color: #868e96; margin-top: 4px; }
      /* Status pills */
      .status-pill { display: inline-block; padding: 4px 12px; border-radius: 12px; font-size: 13px; font-weight: 600; }
      .pill-success { background: #d4edda; color: #155724; }
      .pill-warning { background: #fff3cd; color: #856404; }
      .pill-danger { background: #f8d7da; color: #721c24; }
      /* Recommendation box */
      .rec-box { background: linear-gradient(135deg, #e8f4f8 0%, #ffffff 100%); border: 2px solid #2E86AB; border-radius: 12px; padding: 20px; margin: 15px 0; }
      .rec-box h3 { color: #2E86AB; margin-top: 0; }
      .rec-dose { font-size: 34px; font-weight: 700; color: #2E86AB; }
      /* Loading overlay */
      .shiny-spinner-output-container .load-container { font-size: 14px; }
      /* Button */
      .btn-lg { padding: 14px 24px; font-size: 16px; font-weight: 600; letter-spacing: 0.3px; }
      .btn-primary { background: #2E86AB !important; border-color: #2E86AB !important; }
      .btn-primary:hover { background: #246A8C !important; border-color: #246A8C !important; }
      /* Input labels */
      .control-label, label { font-weight: 500; color: #495057; }
    "))
  ),

  sidebar = sidebar(
    width = 340,
    title = NULL,

    # Patient info section
    div(class = "section-header",
        tags$i(class = "fas fa-user-injured"), "Patient Information"),
    fluidRow(
      column(6, numericInput("age", "Age (yr)", value = 70, min = 18, max = 100)),
      column(6, radioButtons("sex", "Sex", choices = c("M" = 1, "F" = 0),
                              selected = 1, inline = TRUE))
    ),
    fluidRow(
      column(6, numericInput("wt", "Weight (kg)", value = 60, min = 30, max = 150, step = 0.5)),
      column(6, numericInput("ht", "Height (cm)", value = 165, min = 130, max = 210))
    ),
    numericInput("cr", "Serum Cr (mg/dL)", value = 1.0, min = 0.1, max = 10, step = 0.1),

    # Regimen section
    div(class = "section-header",
        tags$i(class = "fas fa-prescription-bottle-medical"), "Current Regimen"),
    numericInput("LD", "Loading dose (mg CBA)", value = 300, min = 100, max = 500, step = 25),
    fluidRow(
      column(6, numericInput("MD", "Maint. dose (mg CBA)", value = 150, min = 25, max = 300, step = 25)),
      column(6, radioButtons("II", "Interval",
                              choices = c("q8h" = 8, "q12h" = 12, "q24h" = 24),
                              selected = 8, inline = TRUE))
    ),
    numericInput("inf_dur", "Infusion duration (hr)", value = 1, min = 0.5, max = 4, step = 0.5),

    # Target
    div(class = "section-header",
        tags$i(class = "fas fa-bullseye"), "Target"),
    numericInput("target_cavg", "Target Cavg total colistin (mg/L)", value = round(TARGET_LOW, 2),
                 min = 1, max = 12, step = 0.25),
    ## 숫자는 상단 상수에서만 온다 — 문구에 하드코딩하면 f_u 가 바뀔 때 화면만 옛 값으로 남는다.
    tags$small(style = "color:#6c757d;",
               HTML(sprintf(paste0(
                 "<b>Aim for %.2f\u2013%.2f mg/L</b> (total drug).<br>",
                 "<b>Lower</b> \u2014 the free-drug target fAUC<sub>24</sub>/MIC \u2265 25 at ",
                 "MIC 1 mg/L, converted using this cohort's measured unbound fraction ",
                 "(f<sub>u,eff</sub> = %.3f).<br>",
                 "<b>Upper</b> \u2014 the 4 mg/L consensus safety limit, converted the same way.<br>",
                 "<span style='color:#a1462b;'>The upper bound is the weaker of the two.</span> ",
                 "The studies behind the 4 mg/L limit never measured unbound fraction, so it could ",
                 "sit anywhere between %.1f and %.1f mg/L \u2014 at the low end, below the efficacy ",
                 "target. Read this window as narrow and approximate, not exact."),
                 TARGET_LOW, TARGET_HIGH, FU_EFF, TARGET_HI_LO, TARGET_HI_HI))),

    # TDM optional
    div(class = "section-header",
        tags$i(class = "fas fa-vial"), "TDM Samples (optional)"),
    materialSwitch("use_tdm", label = "Use TDM for individualization",
                   value = FALSE, status = "primary"),
    conditionalPanel(
      condition = "input.use_tdm == true",
      tags$small("Enter up to 4 TDM observations:", style="color:#6c757d;"),
      tags$hr(style="margin:6px 0;"),
      fluidRow(
        column(4, numericInput("tdm_time1", "Time (hr)", value = 96, min = 0)),
        column(4, selectInput("tdm_cmt1", "Analyte",
                               choices = c("CMS"=1, "ColA"=3, "ColB"=5), selected=3)),
        column(4, numericInput("tdm_dv1", "DV (ng/mL)", value = 8000, min = 0))),
      fluidRow(
        column(4, numericInput("tdm_time2", "Time (hr)", value = 100, min = 0)),
        column(4, selectInput("tdm_cmt2", "Analyte",
                               choices = c("CMS"=1, "ColA"=3, "ColB"=5), selected=3)),
        column(4, numericInput("tdm_dv2", "DV (ng/mL)", value = 12000, min = 0))),
      fluidRow(
        column(4, numericInput("tdm_time3", "Time (hr)", value = 96, min = 0)),
        column(4, selectInput("tdm_cmt3", "Analyte",
                               choices = c("CMS"=1, "ColA"=3, "ColB"=5), selected=5)),
        column(4, numericInput("tdm_dv3", "DV (ng/mL)", value = 800, min = 0))),
      fluidRow(
        column(4, numericInput("tdm_time4", "Time (hr)", value = 100, min = 0)),
        column(4, selectInput("tdm_cmt4", "Analyte",
                               choices = c("CMS"=1, "ColA"=3, "ColB"=5), selected=5)),
        column(4, numericInput("tdm_dv4", "DV (ng/mL)", value = 1100, min = 0)))
    ),

    hr(),
    actionButton("run", tags$span(tags$i(class="fas fa-play"), " Run Analysis"),
                  class = "btn-primary btn-lg w-100")
  ),

  # Main panel
  navset_card_tab(
    id = "main_tabs",
    full_screen = TRUE,

    # ---- Summary tab ----
    nav_panel(
      title = tags$span(tags$i(class="fas fa-chart-line"), " Summary"),

      # Conditional welcome message
      conditionalPanel(
        condition = "input.run == 0",
        div(class = "card",
            div(class = "card-body text-center", style="padding:40px;",
              tags$i(class = "fas fa-hand-point-left",
                     style="font-size:40px; color:#2E86AB; margin-bottom:16px;"),
              tags$h3("Welcome"),
              tags$p("Enter patient information and current regimen on the left, then click ",
                     tags$strong("Run Analysis"), " to begin."),
              tags$p(tags$small("Optional: Enable TDM for Bayesian individualization"))))
      ),
      conditionalPanel(
        condition = "input.run > 0",

        # Patient summary row
        fluidRow(
          column(4, div(class = "metric-box",
            div(class="metric-label", tags$i(class="fas fa-tint"), " CLCR (C-G)"),
            div(class="metric-value", textOutput("metric_clcr", inline=TRUE)),
            div(class="metric-sub", textOutput("metric_clcr_cat", inline=TRUE)))),
          column(4, div(class = "metric-box",
            div(class="metric-label", tags$i(class="fas fa-weight-scale"), " Weight"),
            div(class="metric-value", textOutput("metric_wt", inline=TRUE)),
            div(class="metric-sub", "kg — allometric on CL and V"))),
          column(4, div(class = "metric-box",
            div(class="metric-label", tags$i(class="fas fa-user"), " Patient"),
            div(class="metric-value", textOutput("metric_patient", inline=TRUE)),
            div(class="metric-sub", textOutput("metric_patient_sub", inline=TRUE))))
        ),

        # Current regimen prediction
        br(),
        card(
          card_header(tags$span(tags$i(class="fas fa-flask"), " Current Regimen — Predicted Exposure")),
          card_body(
            withSpinner(uiOutput("current_pred_box"),
                         type = 6, color = "#2E86AB", size = 0.8)
          )
        ),

        # Loading dose Day 0-24h peak warning
        card(
          card_header(tags$span(tags$i(class="fas fa-exclamation-triangle"),
                                 " Loading Dose — Day 0 Peak Check")),
          card_body(
            withSpinner(uiOutput("ld_day1_box"),
                         type = 6, color = "#c62828", size = 0.7)
          )
        ),

        # Recommendation
        withSpinner(uiOutput("recommendation_box"),
                     type = 6, color = "#2E86AB", size = 0.8),

        # PTA table preview
        card(
          card_header(tags$span(tags$i(class="fas fa-table"), " Quick Regimen Comparison")),
          card_body(
            withSpinner(DT::dataTableOutput("dose_options_table_summary"),
                         type = 4, color = "#2E86AB", size = 0.6)
          )
        )
      )
    ),

    # ---- Trajectory tab ----
    nav_panel(
      title = tags$span(tags$i(class="fas fa-wave-square"), " Trajectory"),
      card(
        card_header(tags$span(tags$i(class="fas fa-chart-area"), " 168-hr Concentration-Time Profile")),
        card_body(
          withSpinner(plotOutput("traj_plot", height = "550px"),
                       type = 8, color = "#2E86AB", size = 1)
        )
      ),
      card(
        card_header("Interpretation Guide"),
        card_body(
          tags$ul(
            tags$li(tags$strong("CMS (parent)"), ": Total parent drug"),
            tags$li(tags$strong("Colistin A/B"), ": Active metabolites"),
            tags$li(tags$strong("Total (A+B)"), ": Clinically relevant target; red dashed = 2 mg/L"),
            tags$li("Steady state typically reached at Day 3-5")
          )
        )
      )
    ),

    # ---- Dose Options tab ----
    nav_panel(
      title = tags$span(tags$i(class="fas fa-list"), " Dose Options"),
      card(
        card_header(tags$span(tags$i(class="fas fa-balance-scale"), " All Regimen Options with PTA")),
        card_body(
          withSpinner(DT::dataTableOutput("dose_options_table"),
                       type = 4, color = "#2E86AB", size = 0.8),
          tags$small(tags$i(class="fas fa-info-circle"),
                     " Regimens ranked by PTA (probability of Cavg ≥ target) then median Cavg. Sort by clicking columns.")
        )
      )
    ),

    # ---- TDM Details ----
    nav_panel(
      title = tags$span(tags$i(class="fas fa-vial"), " TDM Details"),
      conditionalPanel(
        condition = "input.use_tdm == true",
        fluidRow(
          column(6, card(
            card_header(tags$span(tags$i(class="fas fa-microscope"), " Estimated Individual ETAs")),
            card_body(
              withSpinner(uiOutput("tdm_etas_ui"),
                           type = 6, color = "#2E86AB", size = 0.8))
          )),
          column(6, card(
            card_header(tags$span(tags$i(class="fas fa-cogs"), " Individual PK Parameters")),
            card_body(
              withSpinner(uiOutput("tdm_params_ui"),
                           type = 6, color = "#2E86AB", size = 0.8))
          ))
        ),
        card(
          card_header(tags$span(tags$i(class="fas fa-info-circle"), " Bayesian Estimation Method")),
          card_body(
            tags$p("MAP (Maximum A Posteriori) Bayesian estimation uses 1000 ETA grid samples from the population OMEGA distribution. Each candidate ETA is evaluated for fit to observed TDM data plus the prior probability. The best-fit ETAs are selected and used to compute individual PK parameters."),
            tags$p(tags$strong("Minimum TDM samples:"), " 2 (ideally ColA + ColB × peak + trough)"),
            tags$p(tags$strong("Best results:"), " Samples at steady state (Day 5+)")
          )
        )
      ),
      conditionalPanel(
        condition = "input.use_tdm == false",
        div(class = "card",
            div(class = "card-body text-center", style="padding:40px;",
              tags$i(class="fas fa-toggle-off", style="font-size:48px; color:#ced4da;"),
              tags$h4("TDM Not Enabled", style="color:#6c757d; margin-top:12px;"),
              tags$p("Toggle ", tags$strong("'Use TDM for individualization'"),
                     " in the sidebar to enable Bayesian personalization."),
              tags$p(tags$small("Without TDM, predictions use population typical values."))))
      )
    ),

    # ---- About ----
    nav_panel(
      title = tags$span(tags$i(class="fas fa-info-circle"), " About"),
      card(
        card_header("Model Information"),
        card_body(
          tags$h5("Six-compartment parent-metabolite model"),
          tags$ul(
            tags$li(tags$strong("Colistimethate (CMS)"), ": two-compartment disposition"),
            tags$li(tags$strong("Colistin A"), ": two-compartment, formed from CMS"),
            tags$li(tags$strong("Colistin B"), ": two-compartment, formed from CMS"),
            tags$li("Each species has its own inter-compartmental clearance. Metabolite volumes are tied to the corresponding CMS volumes (V3 = V5 = V1, V4 = V6 = V2) for identifiability, because no intravenous colistin was given.")
          ),
          tags$h5("Covariates"),
          tags$ul(
            tags$li(tags$strong("Creatinine clearance"), " (Cockcroft-Gault) on all three clearances, as a power function — exponents 0.37, 0.55 and 1.15 for CMS, colistin A and colistin B"),
            tags$li(tags$strong("Body weight"), " by fixed allometric scaling: (WT/70)^0.75 on clearances, (WT/70) on volumes")
          ),
          tags$h5("Population"),
          tags$p("21 adult ICU patients (Korean cohort); 575 observations across the three analytes; creatinine clearance 8.3 - 308.9 mL/min."),
          tags$p("Estimated in NONMEM 7.6.0 with first-order conditional estimation; parameter uncertainty from a 500-replicate nonparametric bootstrap.")
        )
      ),
      card(
        card_header("Limitations"),
        card_body(
          tags$ul(
            tags$li("Developed on N=21; external validation pending"),
            tags$li("Total (bound + unbound) colistin. This cohort's measured unbound fraction is fu,eff = 0.225, well below the 0.5 often assumed — so the exposure shown here is not interchangeable with free-drug targets from cohorts that did not measure binding"),
            tags$li("Intended for adult ICU populations"),
            tags$li("Not a substitute for clinical judgment")
          )
        )
      )
    )
  )
)

# ==================== SERVER ====================
server <- function(input, output, session) {

  derived <- reactive({
    sex_n <- as.numeric(input$sex)
    ## The sex radio maps M -> 1, F -> 0, so the 0.85 Cockcroft-Gault factor belongs on
    ## sex_n == 0. Through Run240 this test was inverted here and in the control stream,
    ## which under-estimated CLCR for every male patient by 15%.
    CLCR <- (140 - input$age)/input$cr * input$wt/72 *
            ifelse(sex_n == 0, 0.85, 1.0)
    CLCR <- min(max(CLCR, 5), 310)   # cohort max under the corrected formula is 308.9
    clcr_cat <- cut(CLCR, c(0,30,60,90,120,310),
                    labels=c("<30 (severe CKD)","30-60","60-90 (normal)",
                             "90-120","ARC (>120)"))
    list(CLCR = round(CLCR, 1),
         CLCR_cat = as.character(clcr_cat))
  })

  etas_est <- reactiveVal(c(0, 0, 0, 0))

  run_map <- function(LD, MD, II, WT, CLCR, observations, N_grid=500) {
    doses <- ev(amt = LD, rate = LD/input$inf_dur, cmt = 1, time = 0) +
             ev(amt = MD, rate = MD/input$inf_dur, ii = II,
                addl = ceiling(168/II) - 1, cmt = 1, time = II)
    etas <- matrix(rnorm(N_grid * 4), N_grid, 4)
    etas <- sweep(etas, 2, OMEGA_SD, "*")
    log_lik <- numeric(N_grid)
    withProgress(message = "Bayesian estimation", value = 0, {
      for (k in 1:N_grid) {
        e <- etas[k, ]
        sim <- mod %>% param(CLCR=CLCR, WT=WT,
                  TVCL_T = pop_par$TVCL_T * exp(e[1]),
                  TVV1   = pop_par$TVV1   * exp(e[2]),
                  TVCLMA = pop_par$TVCLMA * exp(e[3]),
                  TVCLMB = pop_par$TVCLMB * exp(e[4])) %>%
          zero_re() %>% ev(doses) %>%
          mrgsim(end = max(observations$TIME) + 1, delta = 0.5, recsort = 3) %>%
          as.data.frame()
        preds <- sapply(seq_len(nrow(observations)), function(j) {
          r <- observations[j, ]
          pr <- sim %>% filter(abs(time - r$TIME) < 0.5) %>% slice(1)
          if (nrow(pr) == 0) return(NA)
          switch(as.character(r$CMT), "1" = pr$CMS_C, "3" = pr$MA_C, "5" = pr$MB_C, NA)
        })
        sigma_add <- RES_ADD[match(observations$CMT, c(1, 3, 5))]
        sigma_prp <- RES_PRP[match(observations$CMT, c(1, 3, 5))]
        sigma     <- sqrt(sigma_add^2 + (sigma_prp * preds)^2)
        log_lik[k] <- -0.5 * sum(((observations$DV - preds)/sigma)^2 + log(2*pi*sigma^2),
                                  na.rm = TRUE) - 0.5 * sum(e^2)
        if (k %% 50 == 0) incProgress(50/N_grid, detail = sprintf("%d/%d samples", k, N_grid))
      }
    })
    return(etas[which.max(log_lik), ])
  }

  pta_reg <- function(LD, md_test, ii_test, WT, CLCR, etas_center=c(0,0,0,0), use_tdm=FALSE, N_mc=50) {
    set.seed(42)
    etas_mc <- matrix(rnorm(N_mc*4), N_mc, 4)
    etas_mc <- sweep(etas_mc, 2, OMEGA_SD, "*")
    if (use_tdm) etas_mc <- sweep(etas_mc*0.3, 2, etas_center, "+")
    cavgs <- numeric(N_mc)
    cmaxs_cms <- numeric(N_mc)
    for (k in 1:N_mc) {
      ek <- etas_mc[k, ]
      doses_t <- ev(amt = LD, rate = LD/input$inf_dur, cmt=1, time=0) +
                 ev(amt = md_test, rate = md_test/input$inf_dur, ii=ii_test,
                    addl = ceiling(168/ii_test)-1, cmt=1, time=ii_test)
      sim_k <- mod %>% param(CLCR=CLCR, WT=WT,
                  TVCL_T = pop_par$TVCL_T * exp(ek[1]),
                  TVV1   = pop_par$TVV1   * exp(ek[2]),
                  TVCLMA = pop_par$TVCLMA * exp(ek[3]),
                  TVCLMB = pop_par$TVCLMB * exp(ek[4])) %>%
        zero_re() %>% ev(doses_t) %>%
        mrgsim(end=168, delta=2, recsort=3) %>% as.data.frame()
      ss_k <- sim_k %>% filter(time >= 144)
      cavgs[k] <- mean(ss_k$MA_C + ss_k$MB_C, na.rm=TRUE)/1000
      cmaxs_cms[k] <- max(ss_k$CMS_C, na.rm=TRUE)/1000
    }
    list(pta = mean(cavgs >= input$target_cavg, na.rm=TRUE),
         cavg_median = median(cavgs),
         cavg_p25 = quantile(cavgs, 0.25),
         cavg_p75 = quantile(cavgs, 0.75),
         cmax_cms_median = median(cmaxs_cms),
         p_overexposure = mean(cavgs > TARGET_HIGH, na.rm=TRUE))  # P(Cavg > 결합보정 상한)
  }

  results_R <- reactiveValues(
    current = NULL, recommendation = NULL,
    options = NULL, sim = NULL, etas = NULL)

  observeEvent(input$run, {
    d <- derived()
    LD <- input$LD; MD <- input$MD; II <- as.numeric(input$II)

    # TDM if enabled
    if (input$use_tdm) {
      obs <- data.frame(
        TIME = c(input$tdm_time1, input$tdm_time2, input$tdm_time3, input$tdm_time4),
        CMT  = as.numeric(c(input$tdm_cmt1, input$tdm_cmt2, input$tdm_cmt3, input$tdm_cmt4)),
        DV   = c(input$tdm_dv1, input$tdm_dv2, input$tdm_dv3, input$tdm_dv4))
      obs <- obs[!is.na(obs$DV) & obs$DV > 0, ]
      if (nrow(obs) >= 2) {
        etas <- run_map(LD, MD, II, input$wt, d$CLCR, obs, N_grid=1000)
        etas_est(etas)
      } else {
        etas_est(c(0,0,0,0))
      }
    } else {
      etas_est(c(0,0,0,0))
    }
    e <- etas_est()

    # Simulate
    withProgress(message = "Running simulations...", value = 0.3, {
      doses <- ev(amt = LD, rate = LD/input$inf_dur, cmt = 1, time = 0) +
               ev(amt = MD, rate = MD/input$inf_dur, ii = II,
                  addl = ceiling(168/II) - 1, cmt = 1, time = II)
      sim <- mod %>% param(CLCR=d$CLCR, WT=input$wt,
                TVCL_T = pop_par$TVCL_T * exp(e[1]),
                TVV1   = pop_par$TVV1   * exp(e[2]),
                TVCLMA = pop_par$TVCLMA * exp(e[3]),
                TVCLMB = pop_par$TVCLMB * exp(e[4])) %>%
        zero_re() %>% ev(doses) %>%
        mrgsim(end = 168, delta = 1, recsort = 3) %>% as.data.frame()
      results_R$sim <- sim

      ss <- sim %>% filter(time >= 144)
      cavg_curr <- mean(ss$MA_C + ss$MB_C, na.rm = TRUE) / 1000
      cmax_cms  <- max(ss$CMS_C, na.rm = TRUE) / 1000

      # Day 0-24h loading-phase metrics
      day1 <- sim %>% filter(time >= 0, time <= 24)
      day1_cmax_cms <- max(day1$CMS_C, na.rm = TRUE) / 1000
      day1_tot <- (day1$MA_C + day1$MB_C) / 1000
      t_therapeutic <- if (any(day1_tot >= 2)) day1$time[which(day1_tot >= 2)[1]] else NA_real_

      results_R$current <- list(cavg = cavg_curr, cmax_cms = cmax_cms,
                                 day1_cmax_cms = day1_cmax_cms,
                                 t_therapeutic = t_therapeutic,
                                 LD = LD)

      incProgress(0.15, detail = "PTA: current regimen...")
      r_curr <- pta_reg(LD, MD, II, input$wt, d$CLCR, e, input$use_tdm, N_mc=40)
      results_R$current$pta <- r_curr$pta

      # Compute PTA for all common regimens first, THEN pick best recommendation
      options_list <- list(
        c(50, 12), c(75, 12), c(100, 12), c(100, 8),
        c(125, 8), c(150, 8), c(150, 12), c(175, 8), c(200, 8), c(200, 12)
      )
      dose_options <- list()
      total_opts <- length(options_list)
      for (i in seq_along(options_list)) {
        md_opt <- options_list[[i]][1]; ii_opt <- options_list[[i]][2]
        incProgress(0.6/total_opts,
          detail = sprintf("Options %d/%d: %dmg q%dh", i, total_opts, md_opt, ii_opt))
        r <- pta_reg(LD, md_opt, ii_opt, input$wt, d$CLCR, e, input$use_tdm, N_mc=40)
        dose_options[[i]] <- data.frame(
          Dose_mg = md_opt, Interval_hr = ii_opt,
          Daily_mg = round(md_opt * 24 / ii_opt),
          PTA_pct = round(r$pta * 100, 0),
          Cavg_median = round(r$cavg_median, 2),
          Cavg_IQR = sprintf("(%.2f-%.2f)", r$cavg_p25, r$cavg_p75),
          Cmax_CMS_median = round(r$cmax_cms_median, 1),
          P_Cavg_over_pct = round(r$p_overexposure * 100, 0))
      }
      opts_df <- bind_rows(dose_options) %>% arrange(desc(PTA_pct), Cavg_median)
      results_R$options <- opts_df

      # ====================================================
      # Balanced Recommendation (Option A: Tier 1 + Tier 2)
      # ====================================================
      # Tier 1 (must-pass): PTA >= 80%  AND  Cmax_CMS <= 40  AND  P(Cavg > 결합보정 상한) <= 30%
      # Tier 2 (ranking among eligible):
      #   a) Cavg_median in [target, target*3]
      #   b) Lowest overexposure rate
      #   c) Match current interval
      #   d) Lower dose when tied
      # ====================================================
      PTA_THRESHOLD <- 80
      CMAX_CMS_MAX  <- 40
      OVEREXP_MAX   <- 30
      # 창은 결합보정 기준이다(파일 상단 TARGET_* 주석 참조). 사용자의 target_cavg 는
      # 그 안에서의 조준점이며 기본값은 효능 하한 TARGET_LOW (4.63) 이다.
      target_low    <- TARGET_LOW
      target_high   <- TARGET_HIGH

      tier1 <- opts_df %>% filter(PTA_pct >= PTA_THRESHOLD,
                                   Cmax_CMS_median <= CMAX_CMS_MAX,
                                   P_Cavg_over_pct <= OVEREXP_MAX)

      if (nrow(tier1) >= 1) {
        # Rank within Tier 1
        tier1$in_window <- ifelse(tier1$Cavg_median >= target_low &
                                    tier1$Cavg_median <= target_high, 0, 1)
        tier1$interval_match <- ifelse(tier1$Interval_hr == II, 0, 1)
        tier1 <- tier1 %>% arrange(in_window, P_Cavg_over_pct, interval_match, Dose_mg)
        best <- tier1 %>% slice(1)
        rec_status <- "optimal"
        rec_note <- "✅ Meets all safety + efficacy criteria."
      } else {
        # Tier 1 없음: PTA만 달성하는 중에서 Cmax 최소 선택
        tier2_fallback <- opts_df %>% filter(PTA_pct >= PTA_THRESHOLD) %>%
          arrange(Cmax_CMS_median, P_Cavg_over_pct)
        if (nrow(tier2_fallback) >= 1) {
          best <- tier2_fallback %>% slice(1)
          rec_status <- "caution"
          rec_note <- "⚠️ Narrow therapeutic window — close monitoring + TDM recommended."
        } else {
          # PTA 달성도 없음
          best <- opts_df %>% slice(1)  # highest PTA
          rec_status <- "warning"
          rec_note <- "🚨 No regimen meets PTA ≥80%. Individualized TDM essential."
        }
      }

      results_R$recommendation <- list(
        rec_mg = best$Dose_mg, rec_ii = best$Interval_hr,
        cavg_pred = best$Cavg_median, pta = best$PTA_pct / 100,
        cmax_cms = best$Cmax_CMS_median,
        overexp = best$P_Cavg_over_pct,
        status = rec_status, note = rec_note,
        threshold = PTA_THRESHOLD,
        cmax_max = CMAX_CMS_MAX,
        overexp_max = OVEREXP_MAX,
        target_low = target_low, target_high = target_high)
      results_R$etas <- e
      incProgress(0.05, detail = "Done")
    })
  })

  # ===== Outputs =====
  output$metric_clcr <- renderText({ d <- derived(); sprintf("%.0f", d$CLCR) })
  output$metric_clcr_cat <- renderText({ d <- derived(); d$CLCR_cat })
  output$metric_wt  <- renderText({ sprintf("%.0f", input$wt) })
  output$metric_patient <- renderText({
    sprintf("%.0f yr %s", as.numeric(input$age), ifelse(as.numeric(input$sex)==1,"M","F")) })
  output$metric_patient_sub <- renderText({
    sprintf("%.0f kg, Cr %.2f", input$wt, input$cr) })

  output$current_pred_box <- renderUI({
    if (is.null(results_R$current)) return(NULL)
    c <- results_R$current
    pta_class <- if (c$pta >= 0.95) "pill-success"
                 else if (c$pta >= 0.85) "pill-warning" else "pill-danger"
    pta_label <- if (c$pta >= 0.95) "Excellent"
                 else if (c$pta >= 0.85) "Adequate" else "Suboptimal"
    cavg_class <- if (c$cavg >= input$target_cavg && c$cavg <= TARGET_HIGH) "pill-success"
                  else if (c$cavg >= input$target_cavg * 0.5) "pill-warning"
                  else "pill-danger"
    fluidRow(
      column(4, div(class = "metric-box",
        div(class="metric-label", "Cavg Total Colistin"),
        div(class="metric-value", sprintf("%.2f", c$cavg)),
        div(class="metric-sub",
            tags$span(class=paste("status-pill", cavg_class),
              sprintf("Target: %.2f mg/L", input$target_cavg))))),
      column(4, div(class = "metric-box",
        div(class="metric-label", "Cmax CMS"),
        div(class="metric-value", sprintf("%.1f", c$cmax_cms)),
        div(class="metric-sub", "mg/L"))),
      column(4, div(class = "metric-box",
        div(class="metric-label", "PTA (Cavg ≥ target)"),
        div(class="metric-value", sprintf("%.0f%%", c$pta*100)),
        div(class="metric-sub",
            tags$span(class=paste("status-pill", pta_class), pta_label))))
    )
  })

  # Day 0 loading-phase peak check + time-to-therapeutic for user's LD
  output$ld_day1_box <- renderUI({
    if (is.null(results_R$current)) return(NULL)
    c <- results_R$current
    cmax_class <- if (c$day1_cmax_cms >= 80) "pill-danger"
                  else if (c$day1_cmax_cms >= 50) "pill-warning"
                  else "pill-success"
    cmax_label <- if (c$day1_cmax_cms >= 80) "High — consider reduced LD"
                  else if (c$day1_cmax_cms >= 50) "Caution — monitor kidney fn"
                  else "Within typical range"
    t_ther <- c$t_therapeutic
    t_class <- if (is.na(t_ther)) "pill-danger"
               else if (t_ther <= 4) "pill-success"
               else if (t_ther <= 12) "pill-warning" else "pill-danger"
    t_label <- if (is.na(t_ther)) "Not reached in 24 h"
               else if (t_ther <= 4) "Fast (≤ 4 h)"
               else if (t_ther <= 12) "Moderate (4–12 h)"
               else "Slow (> 12 h)"

    tagList(
      tags$p(style = "margin-bottom: 8px; color: #6c757d; font-size: 13px;",
             sprintf("Loading dose: %d mg CBA (1-h IV). Day 0–24 h simulation with individual PK parameters.",
                     c$LD)),
      fluidRow(
        column(6, div(class = "metric-box",
          div(class = "metric-label", "Day 0 peak CMS"),
          div(class = "metric-value", sprintf("%.1f", c$day1_cmax_cms)),
          div(class = "metric-sub",
              tags$span(class = paste("status-pill", cmax_class), cmax_label)))),
        column(6, div(class = "metric-box",
          div(class = "metric-label", "Time to Cavg ≥ 2 mg/L"),
          div(class = "metric-value",
              if (is.na(t_ther)) "— " else sprintf("%.1f h", t_ther)),
          div(class = "metric-sub",
              tags$span(class = paste("status-pill", t_class), t_label))))
      ),
      tags$p(style = "margin-top: 8px; color: #6c757d; font-size: 12px;",
             "Reference: CMS Cmax > 50 mg/L considered cautionary (extrapolated from safety PK literature); > 80 mg/L warrants reducing LD or extending infusion duration. The TDM tool does not automatically adjust the loading dose — TDM sampling occurs after the loading phase. Clinical judgement advised.")
    )
  })

  output$recommendation_box <- renderUI({
    if (is.null(results_R$recommendation)) return(NULL)
    r <- results_R$recommendation
    cur_MD <- as.numeric(input$MD); cur_II <- as.numeric(input$II)

    direction_icon <- if (r$rec_mg > cur_MD) "fa-arrow-up"
                      else if (r$rec_mg < cur_MD) "fa-arrow-down" else "fa-equals"
    direction_color <- if (r$rec_mg > cur_MD) "#dc3545"
                       else if (r$rec_mg < cur_MD) "#28a745" else "#6c757d"
    direction_text <- if (r$rec_mg > cur_MD) "INCREASE"
                      else if (r$rec_mg < cur_MD) "DECREASE" else "KEEP"
    ii_changed_text <- if (r$rec_ii != cur_II) sprintf(" + interval → q%.0fh", r$rec_ii) else ""

    # Traffic light for each criterion
    pta_light <- if (r$pta*100 >= 90) "🟢"
                 else if (r$pta*100 >= r$threshold) "🟡" else "🔴"
    cavg_light <- if (r$cavg_pred >= r$target_low & r$cavg_pred <= r$target_high) "🟢"
                  else if (r$cavg_pred < r$target_low) "🔴"
                  else if (r$cavg_pred <= 8) "🟡" else "🔴"
    cmax_light <- if (r$cmax_cms <= r$cmax_max) "🟢"
                  else if (r$cmax_cms <= r$cmax_max*1.3) "🟡" else "🔴"
    overexp_light <- if (r$overexp <= r$overexp_max) "🟢"
                     else if (r$overexp <= r$overexp_max*1.5) "🟡" else "🔴"

    status_banner <- switch(r$status,
      "optimal" = tags$div(style="background:#d4edda; color:#155724; padding:8px 14px; border-radius:6px; font-weight:600; margin-bottom:10px;",
                            tags$i(class="fas fa-check-circle"),
                            " OPTIMAL: Meets all efficacy + safety criteria"),
      "caution" = tags$div(style="background:#fff3cd; color:#856404; padding:8px 14px; border-radius:6px; font-weight:600; margin-bottom:10px;",
                            tags$i(class="fas fa-exclamation-triangle"),
                            " CAUTION: Narrow therapeutic window — TDM recommended"),
      "warning" = tags$div(style="background:#f8d7da; color:#721c24; padding:8px 14px; border-radius:6px; font-weight:600; margin-bottom:10px;",
                            tags$i(class="fas fa-times-circle"),
                            " WARNING: No ideal regimen — individualized TDM essential"))

    div(class = "rec-box",
      status_banner,
      tags$h3(tags$i(class="fas fa-star"), " Recommended Dose"),
      fluidRow(
        column(5,
          div(class = "rec-dose",
              sprintf("%.0f mg q%.0fh", r$rec_mg, r$rec_ii)),
          tags$div(style=sprintf("color:%s; font-weight:600; margin-top:8px;", direction_color),
            tags$i(class=paste("fas", direction_icon)),
            sprintf(" %s from %.0f mg q%.0fh%s",
                    direction_text, cur_MD, cur_II, ii_changed_text))),
        column(7,
          # Traffic lights
          tags$table(style="width:100%; border-collapse:collapse;",
            tags$tr(
              tags$td(style="padding:4px 8px;", pta_light, tags$strong(" PTA:")),
              tags$td(style="padding:4px 8px;", sprintf("%.0f%%", r$pta*100)),
              tags$td(style="padding:4px 8px; color:#6c757d; font-size:12px;",
                      sprintf("≥%d%% target", r$threshold))),
            tags$tr(
              tags$td(style="padding:4px 8px;", cavg_light, tags$strong(" Cavg:")),
              tags$td(style="padding:4px 8px;", sprintf("%.2f mg/L", r$cavg_pred)),
              tags$td(style="padding:4px 8px; color:#6c757d; font-size:12px;",
                      sprintf("window %.1f–%.1f", r$target_low, r$target_high))),
            tags$tr(
              tags$td(style="padding:4px 8px;", cmax_light, tags$strong(" Cmax CMS:")),
              tags$td(style="padding:4px 8px;", sprintf("%.1f mg/L", r$cmax_cms)),
              tags$td(style="padding:4px 8px; color:#6c757d; font-size:12px;",
                      sprintf("≤ %.0f mg/L", r$cmax_max))),
            tags$tr(
              tags$td(style="padding:4px 8px;", overexp_light, tags$strong(" Overexp:")),
              tags$td(style="padding:4px 8px;", sprintf("%.0f%%", r$overexp)),
              tags$td(style="padding:4px 8px; color:#6c757d; font-size:12px;",
                      sprintf("\u2264 %d%% (Cavg > %.2f)", r$overexp_max, TARGET_HIGH)))
          )
        )
      ),
      tags$hr(),
      # Note + criteria explanation
      tags$div(
        tags$p(style="margin-bottom:6px;",
               tags$strong("Status: "), r$note),
        tags$details(
          tags$summary(style="cursor:pointer; color:#2E86AB; font-weight:600;",
                       tags$i(class="fas fa-list-check"), " Decision criteria (click to expand)"),
          tags$div(style="margin-top:10px; padding:10px; background:#f8f9fa; border-radius:6px; font-size:13px;",
            tags$p(tags$strong("Tier 1 — Must-pass (all required):")),
            tags$ul(
              tags$li(sprintf("Efficacy: PTA ≥ %d%% for Cavg ≥ %.1f mg/L", r$threshold, input$target_cavg)),
              tags$li(sprintf("Safety (parent): Median Cmax CMS ≤ %.0f mg/L", r$cmax_max)),
              tags$li(sprintf("Safety (metabolite): P(Cavg > %.2f mg/L) \u2264 %d%% \u2014 binding-adjusted upper bound, f_u,ref = 0.31 (sensitivity range %.1f\u2013%.1f)",
                              TARGET_HIGH, r$overexp_max, TARGET_HI_LO, TARGET_HI_HI))),
            tags$p(tags$strong("Tier 2 — Ranking among Tier 1 passed:")),
            tags$ol(
              tags$li(sprintf("Cavg in therapeutic window [%.1f, %.1f] mg/L", r$target_low, r$target_high)),
              tags$li(sprintf("Lowest overexposure risk (P(Cavg > %.2f mg/L))", TARGET_HIGH)),
              tags$li("Match current interval (clinical convenience)"),
              tags$li("Lower dose when tied (cost/toxicity)")),
            tags$p(style="font-size:11px; color:#6c757d; margin-top:8px;",
                   tags$i(class="fas fa-info-circle"),
                   sprintf(" Therapeutic window = [target, target×3] = [%.1f, %.1f]. Loading dose %.0f mg unchanged.",
                           r$target_low, r$target_high, as.numeric(input$LD)))
          )
        )
      )
    )
  })

  output$traj_plot <- renderPlot({
    if (is.null(results_R$sim)) return(NULL)
    sim_long <- results_R$sim %>%
      select(time, CMS_C, MA_C, MB_C) %>%
      mutate(Total = MA_C + MB_C) %>%
      pivot_longer(c(CMS_C, MA_C, MB_C, Total), names_to="Analyte", values_to="Conc_ngmL") %>%
      mutate(Conc_mgL = Conc_ngmL/1000,
             Analyte = factor(Analyte, levels=c("CMS_C","MA_C","MB_C","Total"),
               labels=c("CMS (parent)","Colistin A","Colistin B","Total Colistin (A+B)")))
    ggplot(sim_long, aes(time, Conc_mgL, color=Analyte)) +
      geom_line(linewidth=1.1) +
      geom_hline(yintercept=input$target_cavg, linetype="dashed", color="#e63946",
                 linewidth=0.7) +
      annotate("text", x=155, y=input$target_cavg*1.3,
               label=sprintf("Target %.0f mg/L", input$target_cavg),
               color="#e63946", size=3.5) +
      scale_y_log10() +
      scale_color_manual(values=c("CMS (parent)"="#2E86AB",
                                   "Colistin A"="#A23B72",
                                   "Colistin B"="#F18F01",
                                   "Total Colistin (A+B)"="#C73E1D")) +
      facet_wrap(~Analyte, scales="free_y") +
      labs(title=sprintf("Predicted Concentration-Time (LD %.0f + %.0f q%.0fh)",
                          as.numeric(input$LD), as.numeric(input$MD), as.numeric(input$II)),
           subtitle=ifelse(input$use_tdm,
             "▸ With Bayesian-estimated individual ETAs",
             "▸ Population typical prediction (no TDM)"),
           x="Time (hr)", y="Concentration (mg/L, log)") +
      theme_minimal(base_size=13) +
      theme(legend.position="none",
            plot.title=element_text(face="bold", color="#2E86AB"),
            strip.background=element_rect(fill="#f8f9fa", color=NA),
            strip.text=element_text(face="bold"))
  })

  output$dose_options_table <- DT::renderDataTable({
    if (is.null(results_R$options)) return(NULL)
    DT::datatable(results_R$options,
      options=list(pageLength=14, dom='ft',
                    columnDefs=list(list(className='dt-center', targets='_all'))),
      rownames=FALSE,
      colnames=c("Dose (mg)", "Interval (hr)", "Daily (mg)",
                 "PTA (%)", "Cavg median", "Cavg IQR")) %>%
      DT::formatStyle("PTA_pct",
        background = DT::styleColorBar(c(0, 100), "#2E86AB"),
        backgroundSize = "100% 70%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "right center",
        color = 'white', fontWeight = 'bold') %>%
      DT::formatStyle("Cavg_median",
        backgroundColor = DT::styleInterval(c(2, 4, 8),
          c("#f8d7da", "#d4edda", "#fff3cd", "#f8d7da")))
  })

  output$dose_options_table_summary <- DT::renderDataTable({
    if (is.null(results_R$options)) return(NULL)
    DT::datatable(head(results_R$options, 5),
      options=list(dom='t', ordering=FALSE),
      rownames=FALSE,
      colnames=c("Dose (mg)", "Interval", "Daily", "PTA (%)", "Cavg", "IQR"),
      caption="Top 5 options (by PTA)") %>%
      DT::formatStyle("PTA_pct",
        background = DT::styleColorBar(c(0, 100), "#2E86AB"),
        backgroundSize = "100% 70%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "right center",
        color='white', fontWeight='bold')
  })

  output$tdm_etas_ui <- renderUI({
    if (is.null(results_R$etas)) return(NULL)
    e <- results_R$etas
    div(
      div(class="metric-box",
          div(class="metric-label", "ETA CL"), div(class="metric-value", sprintf("%.3f", e[1]),
              style="font-size:20px;")),
      div(class="metric-box",
          div(class="metric-label", "ETA V1"), div(class="metric-value", sprintf("%.3f", e[2]),
              style="font-size:20px;")),
      div(class="metric-box",
          div(class="metric-label", "ETA CLMA"), div(class="metric-value", sprintf("%.3f", e[3]),
              style="font-size:20px;")),
      div(class="metric-box",
          div(class="metric-label", "ETA CLMB"), div(class="metric-value", sprintf("%.3f", e[4]),
              style="font-size:20px;")))
  })

  output$tdm_params_ui <- renderUI({
    if (is.null(results_R$etas)) return(NULL)
    d <- derived(); e <- results_R$etas
    ## Must mirror $MAIN of the loaded model exactly — this panel is what a clinician reads
    ## off the screen, so a stale exponent here is a silently wrong reported clearance.
    WTR  <- input$wt/70
    CL   <- pop_par$TVCL_T  * (d$CLCR/70)^EXP_CL * WTR^EXP_WT_C * exp(e[1])
    V1   <- pop_par$TVV1    * WTR^EXP_WT_V * exp(e[2])
    CLMA <- pop_par$TVCLMA  * (d$CLCR/70)^EXP_MA * WTR^EXP_WT_C * exp(e[3])
    CLMB <- pop_par$TVCLMB  * (d$CLCR/70)^EXP_MB * WTR^EXP_WT_C * exp(e[4])
    div(
      div(class="metric-box",
          div(class="metric-label", "CL_CMS"), div(class="metric-value", sprintf("%.3f", CL),
              style="font-size:20px;"), div(class="metric-sub", "L/hr")),
      div(class="metric-box",
          div(class="metric-label", "V1"), div(class="metric-value", sprintf("%.2f", V1),
              style="font-size:20px;"), div(class="metric-sub", "L")),
      div(class="metric-box",
          div(class="metric-label", "CLMA"), div(class="metric-value", sprintf("%.3f", CLMA),
              style="font-size:20px;"), div(class="metric-sub", "L/hr")),
      div(class="metric-box",
          div(class="metric-label", "CLMB"), div(class="metric-value", sprintf("%.3f", CLMB),
              style="font-size:20px;"), div(class="metric-sub", "L/hr")))
  })
}

shinyApp(ui = ui, server = server)
