# CV 修改指南 / CV Modification Guide

## 问题：CV 那部分在哪里修改来着？

您的学术网站使用 al-folio 主题，有**两种**方式来修改 CV 内容。系统会按优先级选择数据源：

## 🎯 当前系统工作方式

1. **优先级 1**：如果 `assets/json/resume.json` 存在 → 使用 JSON Resume 格式
2. **优先级 2**：如果 JSON 文件不存在 → 回退到 `_data/cv.yml` (YAML 格式)

**当前状态**：您的网站同时拥有两个文件，所以正在使用 **JSON 格式**。

## 方式 1：JSON Resume 格式 ⭐ (当前使用)

### 文件位置：
- **数据文件**：`assets/json/resume.json` 
- **配置文件**：`_config.yml` (第628-640行)

### 特点：
- ✅ 遵循 [JSON Resume 国际标准](https://jsonresume.org/) 
- ✅ 结构化、专业化
- ✅ 支持导出为其他格式
- ✅ 当前正在使用

### 当前数据结构：
```json
{
  "basics": {
    "name": "Jian Wang", 
    "label": "PhD Student in Theoretical Microbial Ecology",
    "email": "jian.wang1@kuleuven.be",
    "location": { "city": "Ghent", "countryCode": "BE" }
  },
  "education": [
    {
      "institution": "KU Leuven",
      "area": "Theoretical Microbial Ecology", 
      "studyType": "PhD",
      "startDate": "2022-08",
      "endDate": "2026-08"
    }
  ],
  "work": [...],
  "awards": [...],
  "skills": [...],
  "projects": [...]
}
```

### 如何修改：
1. 直接编辑 `assets/json/resume.json`
2. 按照 JSON Resume 标准修改内容
3. 推送到 GitHub 或本地重新构建

## 方式 2：YAML 格式 (备用方案)

### 文件位置：
- **数据文件**：`_data/cv.yml`

### 特点：
- ✅ 更易读易编辑
- ✅ 支持自定义显示类型
- ❌ 仅在删除 JSON 文件后生效

### 当前数据结构：
```yaml
- title: General Information
  type: map
  contents:
    - name: Full Name
      value: Jian Wang
    - name: Date of Birth
      value: 04th May 1994

- title: Education  
  type: time_table
  contents:
    - title: PhD-ing in Theoretical Microbial Ecology
      institution: KU Leuven, Ghent, Belgium
      year: August 2022 - August 2026
      description:
        - Research focused on developing mathematical tools...
```

### 支持的显示类型：
- `map`: 键值对列表
- `time_table`: 时间线表格  
- `list`: 简单列表
- `nested_list`: 嵌套列表
- `list_groups`: 分组列表

### 如何切换到 YAML 方式：
1. 删除或重命名 `assets/json/resume.json`
2. 编辑 `_data/cv.yml` 
3. 重新构建网站

## 📁 其他相关文件

### 页面配置
- **文件**：`_pages/cv.md`
- **作用**：CV 页面基本设置
- **可修改**：页面标题、描述、PDF 下载链接

### 样式文件  
- **文件**：`_sass/_cv.scss`
- **作用**：CV 页面外观样式
- **可修改**：颜色、字体、布局等

### 模板文件
- **文件**：`_layouts/cv.liquid` 
- **作用**：CV 页面 HTML 结构
- **高级用户**：可修改页面布局逻辑

### 配置文件
- **文件**：`_config.yml` (第627-640行)
- **作用**：JSON Resume 数据获取和显示配置

## 🚀 快速修改步骤

### 如果使用 JSON 格式 (推荐)：
1. 编辑 `assets/json/resume.json`
2. 修改对应的 JSON 字段
3. 保存并推送到 GitHub

### 如果切换到 YAML 格式：
1. 删除 `assets/json/resume.json`
2. 编辑 `_data/cv.yml`  
3. 保存并推送到 GitHub

## 🔄 重新构建网站

### 本地开发：
```bash
bundle exec jekyll serve --host 0.0.0.0 --port 4000
```

### GitHub Pages：
推送到 GitHub 后自动重新构建

## 💡 建议

1. **继续使用 JSON 格式**：如果您喜欢标准化的专业格式
2. **切换到 YAML 格式**：如果您更喜欢易读的格式和自定义显示类型
3. **保持数据同步**：如果保留两个文件，确保内容一致

---

*💡 提示：JSON 格式遵循国际标准，便于与其他平台集成；YAML 格式提供更多显示选项，便于本地编辑。*