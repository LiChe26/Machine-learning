# Depression predictor — GitHub + Posit Connect Cloud

这是部署准备包，并非已经上线的应用。

## 已做的部署调整

- 将入口文件命名为 app.R。
- 删除未使用的训练/分析包，显式加载 mlr3learners 和 xgboost。
- 保留原数据、参数、界面及预测逻辑。
- setup.R 在你的 R 环境生成实际安装版本的 manifest.json；不伪造依赖版本。

## 操作步骤（Mac）

1. 下载并解压 depression_github_deploy.zip，将 depression-predictor 文件夹放在 Downloads 下。
2. 打开 RStudio，在 Console 运行：

```r
setwd(path.expand("~/Downloads/depression-predictor"))
source("setup.R")
shiny::runApp(".")
```

3. 页面打开后点 RUN，检查概率图和 Local SHAP。若报错，保存完整报错，先解决再发布。
4. 浏览器打开 https://github.com/new 。Repository name 填 depression-predictor。免费 Connect Cloud 的 GitHub 部署需要 Public 仓库；确认嵌入的数据允许公开。创建仓库。
5. 在仓库页面用 uploading an existing file 或 Add file → Upload files，将 app.R、manifest.json、README.md 拖入仓库根目录。不要上传 zip 本身。点击 Commit changes。
6. 打开 https://connect.posit.cloud/ 并登录。点击 Publish，按提示连接 GitHub，仅授权相应仓库。选择 Shiny、depression-predictor 仓库、main 分支、app.R 主文件，点击 Publish。
7. 部署成功后，复制平台显示的应用链接。最终访问地址由平台分配，无法在部署前确定。

## 必须核对的计算问题

原代码训练数据中的 Age、PIR、ALB、HR 为原始数值，但新输入在预测时被转换为 z 分数，训练和预测的尺度不一致。本包未修改这一计算逻辑。正式作为论文工具或供他人使用前，须根据研究中最终模型的预处理流程纠正并核对结果。

当前每次点击 RUN 都重新训练模型。正式应用建议加载固定的最终模型，核对与论文一致；同时避免重复训练导致等待和结果波动。

## 验证范围

已检查源文件、依赖调用和打包内容。当前制作环境没有 R，未执行 R 安装、应用运行或云端部署。本包需要在你的 RStudio 中完成上述运行检查。

## 官方操作依据

- https://docs.posit.co/connect-cloud/how-to/r/shiny-r.html
- https://docs.github.com/en/pages/getting-started-with-github-pages/what-is-github-pages
