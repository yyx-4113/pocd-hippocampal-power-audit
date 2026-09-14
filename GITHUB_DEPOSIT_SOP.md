# 剩余操作手册（纯中文）

更新时间：2026-09-14 16:06

> **阅读约定**：本手册正文全中文。GitHub 与 Zenodo 网页默认是英文界面，
> 因此每步末尾单独附一张「界面英文对照表」，只在那里出现英文单词。
>
> **更省事的做法**：用 Edge 或 Chrome 打开这些网页后，在页面空白处
> **右键 → 翻译成中文**，界面会变成中文，下面的步骤可以完全按中文照做。

---

## 已经做好的（不用你动手）

| 项目 | 状态 |
|---|---|
| 复现包组装 | 已在本地完成：203 个文件，77.7 MB，含全部脚本、中间表、图与 GEO 缓存 |
| `README.md` / `CITATION.cff` / `LICENSE`（MIT）/ `.gitignore` | 已写好 |
| `.github/workflows/release.yml` | 已放好：打 tag 会自动建 Release |
| `author_verification_statement.md` | 已写好（核验内容与日期见文内） |
| `step5d_composition_animallevel.py` | 新写：把稿件中心表从"无脚本可复现"补成 `python 一条命令` 可复现，已核对与归档表**逐位一致** |
| 用户名纠正 | 账号确认是 **yyx-4113**（旧的 `yongxinyang` 是早年课程作业账号，勿再用） |

下面 3 件事**只能你本人做**（需要账号权限，本机够不着网页端）。

---

# 第一步：新建仓库（约 1 分钟，必须最先做）

推送需要仓库先存在；SSH 密钥能让本机往**已存在**的仓库推代码，但不能替你建仓库。

**1. 打开新建页面**

```
https://github.com/new
```

**2. 填三处**

| 字段 | 填什么 |
|---|---|
| Repository name | `pocd-hippocampal-power-audit`（必须完全一致，含连字符） |
| Description | 见下方那一整行，复制粘贴 |
| 可见性 | 选 **Public** |

描述粘这一整行（英文，面向国外读者，勿改中文）：

```
Reproducibility bundle for an animal-level power audit of a 3-vs-3 hippocampal single-cell RNA-seq study in postoperative neurocognitive dysfunction (unit of inference, pseudobulk state tests, ambient-RNA attribution, required-n calibration).
```

**3. 右上角那些勾全部不要勾**

README、.gitignore、License 三个都留空 —— 我们已经准备好了，
如果 GitHub 替你建了 README，推送时会冲突，还得再多一步。

**4. 点最下面绿色的创建按钮。**

**5. 建好后跟我说一声**，我直接把 203 个文件推上去、打 v1.0.0 标签、
发布页会自动建好。

### 出问题了怎么办

| 现象 | 处理 |
|---|---|
| 提示名字已被占用 | 说明你账号下已有同名仓库，把名字后面加 `-v1` 告诉我 |
| 建完发现是 Private | 到 `设置 → 通用 → 最下方 Danger Zone → 改可见性`，改成 Public |
| 不小心勾了 README | 告诉我，我推送前先合并一次历史，不影响结果 |

### 界面英文对照表

| 中文意思 | 网页上的英文原词 |
|---|---|
| 新建仓库 | New repository |
| 仓库名 | Repository name |
| 描述 | Description |
| 公开 | Public |
| 私有 | Private |
| 创建仓库 | Create repository |

---

# 第二步：在 Zenodo 上获取永久编号（约 3 分钟，可在推送后做）

**1. 登录 Zenodo**

打开 `https://zenodo.org`，点右上角登录，选**带 GitHub 图标的那一项**。
（如果本机打不开 Zenodo，见文末备选方案。）

**2. 打开仓库开关**

```
https://zenodo.org/account/settings/github/
```

在列表里找到 `yyx-4113/pocd-hippocampal-power-audit`，把这一行**最右边的开关**点开。

**3. 等它归档**，通常几分钟。想催一下：回到仓库的 Release 页 →
编辑 → 什么都不改 → 更新发布。

**4. 复制编号**

刷新上面那个设置页，点进这条记录，复制 `10.5281/zenodo.XXXXXXX` 这一串**发我**。

> 坑：详情页上如果有个开关叫「引用所有版本」，**保持关闭**，否则编号会随版本变动。

**5. 我拿到编号后一次性做完**：回填 `README.md` 与
`author_verification_statement.md` 的 DOI 占位、回填稿件 Data availability、
重建 Word 文档、同步到仓库。

### 界面英文对照表

| 中文意思 | 网页上的英文原词 |
|---|---|
| 登录 | Log in |
| 用 GitHub 登录 | Log in with GitHub |
| 授权 | Authorize zenodo |
| 立即同步 | sync now |
| 引用所有版本 | Cite all versions? |
| 更新发布 | Update release |

---

# 第三步：把默认分支确认为 main（约 30 秒，建议做）

新仓库默认分支一般就是 `main`，但 GitHub 有时仍给 `master`。
打开：

```
https://github.com/yyx-4113/pocd-hippocampal-power-audit/settings/branches
```

确认「默认分支」是 `main`；不是的话点右侧双向箭头切换成 `main` 再保存。

### 界面英文对照表

| 中文意思 | 网页上的英文原词 |
|---|---|
| 默认分支 | Default branch |
| 双向箭头 | ⇄ 图标 |
| 更新 | Update |

---

# 附：本机环境备忘（供以后复用）

- 本机 `github.com` 网页可达，`api.github.com` 可达，`zenodo.org` 测试不通；
  推送走 SSH（`git@github.com`），密钥为 `~/.ssh/id_ed25519`，身份 `yyx-4113`
- `git` 需用绝对路径调用（`/c/Program Files/Git/bin/git.exe`），直接敲 `git` 在本沙箱会挂起
- `gh` 已安装但未登录；本流程不依赖它登录（建 Release 由仓库内的 workflow 完成）
- 未能自动化：**新建仓库**、改仓库描述、改默认分支 —— 这三项只能本人在网页端做

## 备选方案：Zenodo 打不开时改用科学数据银行

1. 打开 `https://www.scidb.cn` → 注册登录
2. 新建数据集 → 填标题、作者（Yang Yongxin）、摘要
3. 上传打包好的源码，浏览器直接下载：

   ```
   https://github.com/yyx-4113/pocd-hippocampal-power-audit/archive/refs/tags/v1.0.0.tar.gz
   ```

4. 提交 → 审核通过后拿到编号，发我即可
