## deploy.R — TDM 앱을 shinyapps.io 에 올린다.
##
## ⚠️ 먼저 읽을 것
##   이 앱의 **치료 창이 바뀌었고, 그래서 권고 용량이 달라진다.** 이전 판은 폐기된
##   총약물 2–8 mg/L 창을 라벨뿐 아니라 용량 선택 로직에도 쓰고 있었고, 지금은
##   결합보정 창(효능 4.63 / 안전 5.51 mg/L 총약물)을 쓴다. 같은 환자에게
##   150 mg q8h 를 통과시키던 것이 이제 100 mg q12h 로 감량을 권고한다.
##   **공개 배포는 저자 합의 뒤에 한다.**
##
## 사용법
##   1) 자격증명을 한 번만 설정한다 (shinyapps.io > Account > Tokens 에서 복사).
##      비밀번호·토큰은 이 파일에 적지 말 것 — 아래는 대화형으로 붙여넣는 방식이다.
##        rsconnect::setAccountInfo(name="<계정>", token="<토큰>", secret="<시크릿>")
##   2) 이 폴더에서 실행:
##        Rscript deploy.R
##
## 앱 이름
##   기본은 새 이름 'colistin-tdm' 이다. 2023년에 올라간 'colistin_new' 는
##   **다른(구) 모형**이므로 덮어쓰지 않는다. 구버전을 내리려면 아래 §3 참고.

library(rsconnect)

APP_NAME <- "colistin-tdm"          # 필요하면 바꾼다
FILES    <- c("app.R", "colistin_run330_mrgsolve.rds")

stopifnot(all(file.exists(FILES)))
if (nrow(accounts()) == 0)
  stop("shinyapps.io 계정이 설정되지 않았다. 위 1) 단계를 먼저 실행할 것.")

cat("배포 대상:", APP_NAME, "\n")
cat("번들 파일:", paste(FILES, collapse = ", "), "\n")
cat("계정      :", accounts()$name[1], "\n\n")

deployApp(appDir   = ".",
          appFiles = FILES,
          appName  = APP_NAME,
          appTitle = "Colistin Bayesian TDM (Run330)",
          forceUpdate = TRUE)

## §3. 2023년 구버전('colistin_new')을 내리려면 — 잘못된 모형을 보여주고 있으므로 권장
##     rsconnect::terminateApp("colistin_new")
