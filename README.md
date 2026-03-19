# Алгоритмы искусственного интеллекта в биологии и медицине

Репозиторий курса **ДПП ПП «Основы биоинформатики и компьютерной биологии»**  
ДВ3: Алгоритмы ИИ в биологии и медицине

---

## 📚 Содержание курса

### Тема 1. Основы построения моделей машинного обучения для биомедицинских данных

**Лекция:** основы классификации и кластеризации, типы биомедицинских данных, препроцессинг, метрики качества (Accuracy, AUC-ROC, F1), переобучение и недообучение.
Видео: 
https://rutube.ru/video/private/71da492db9a390317c68dc52884c878a/?p=NN1VcvV_IdGjMb8JJW10Sg
https://rutube.ru/video/private/945c5a750834ee5302833347ee5257c0/?p=y65fSTpIfSc3OO9WMlVolA

**Практика:**
- `Практика_1-1.R` — логистическая регрессия на Breast Cancer Wisconsin
  Видео:
  https://rutube.ru/video/private/f6653f2631bd4aaac5664839edcd2e04/?p=jubjdr9bAghKk50PZlc5Zg
  https://rutube.ru/video/private/1287ffaeeb608aee52a57955eb9125fe/?p=6rPH1kyctIHNjvNH4Akpnw
- `Практика_1-2.R` — сравнение классификаторов (Logistic Regression, Decision Tree, Random Forest)
  Видео:
  https://rutube.ru/video/private/ec7c811716bd4926f22315898afaa7fe/?p=2lDzrtAPp1NCVQwsKOD_cQ
  https://rutube.ru/video/private/7e5ac56857ecddf52d461d858b632ab4/?p=qL0QXkd3uNnf4rrQ8wQc2g
- `Практика_1-3.R` — кластеризация, оценка числа кластеров
  Видео:
  https://rutube.ru/video/private/359bd930c0e2f9bedbf6915e4711a655/?p=oSEyoTd3P59yse2iwnwYdg
  https://rutube.ru/video/private/e680ad998c1944b0b112576c2b42f174/?p=T_24WAK2cX2fRIPFC4YR7Q

**Задания:** 1–2 (файл с заданиями в корне репозитория)

---

### Тема 2. Продвинутое моделирование и оптимизация

**Лекция:** Feature Importance, SHAP/LIME, работа с несбалансированными данными (SMOTE), гиперпараметрическая оптимизация (Grid Search, Random Search, Байесовская оптимизация), ансамблевые методы.
Видео:
https://rutube.ru/video/private/dad1884fb675b6593909c46640909b70/?p=qJZokxsam1oLrrcdplTv9A
https://rutube.ru/video/private/2255d448364e404bbf70ce339b8f43bd/?p=yricLqVRN3Ku8jaY2iME8Q

**Практика:**
- `Практика 2-1.R` — отбор признаков (Feature Selection), топ-10 генов
  Видео:
  https://rutube.ru/video/private/1dc4666b5a461db8b38d786e0a15198a/?p=qwwdvi8G7ktk8gREZweNNw
  https://rutube.ru/video/private/a5149e0089c1145bc840f0ea57b6220d/?p=hfRFVXyjvdR8YRQ8S35PLQ
- `Практика 2-2.R` — Grid Search оптимизация Random Forest
  Видео:
  https://rutube.ru/video/private/7c4b5014fba4e3d29debcba4c1b4a841/?p=MdvM22fHhsc8UJS2VHlPPQ
  https://rutube.ru/video/private/beae70bc43828a0626870c8225c09505/?p=ELIl8ta8R7Yl1638PwbI_Q
- `Практика 2-3.R` — работа с несбалансированными данными, SMOTE

**Задания:** 3–4

---

### Тема 3. Построение и интерпретация моделей выживаемости *(в разработке)*

**Лекция:** кривые Каплана–Мейера, лог-ранговый тест, модель пропорциональных рисков Кокса, Hazard Ratio, C-index, конкурирующие риски, time-dependent ковариаты, LASSO/Ridge Cox.
Видео:
https://rutube.ru/video/private/8b1ed8e2a3cdb399ed81018db904d73b/?p=v5SxAIdYyL1n70ILppQ4VA
https://rutube.ru/video/private/f1c0a1f73128de1c5c52b1d281f29060/?p=q0J0mgK4mWj3kF6FX1Bt_Q
https://rutube.ru/video/private/bf23c8680008ec83f13530a9979f81c7/?p=JDnAZjM_qHY7anUiqd8TMg

**Практика:**
- `Практика 3-1.r` — кривые Каплана–Мейера, лог-ранговый тест
- `Практика 3-2.r` — одномерная и многомерная модель Кокса, Hazard Ratio
- `Практика 3-3.r` — продвинутые модели выживаемости (конкурирующие риски, time-dependent ковариаты, LASSO Cox)

**Задания:** 5–6

---

### Тема 4. Сетевой анализ и комплексная интерпретация биологических систем *(в разработке)*

**Лекция:** WGCNA, теория графов, обнаружение модулей генов, хаб-гены (kME), функциональное обогащение (GO, KEGG), интеграция с ML (WGCNA + Cox / Random Forest).
Видео:
https://rutube.ru/video/private/6aabd284f97168ecd9d65d8b5f4510c5/?p=OoawCTr4UeS64DnCdy1EPA
https://rutube.ru/video/private/0a4549cbf2eb0a8149e300ec68eedb8c/?p=5EyjMAVtqUatJEk2vaBVDA
https://rutube.ru/video/private/ecb65eee35adbadc019fa082d46309b1/?p=k_8Qg48B2AOt3TF4_uq5Jg

**Практика:** в разработке

**Задания:** 7–8

---

### Тема 5. Нейронные сети в биомедицине

**Лекция:** Perceptron, MLP, CNN (гистология, МРТ), RNN/LSTM (ЭКГ, ЭЭГ), Transfer Learning, GAN (аугментация медицинских изображений), GNN (предсказание взаимодействий лекарство–мишень), AlphaFold, интерпретируемость (GradCAM, SHAP).
Видео:
https://rutube.ru/video/private/9043de7872b1ca9b0fe6e7d0c2b0c136/?p=B7x8UFTCe7X0cKQKiB883A
https://rutube.ru/video/private/b287727175d03716e8f54b7c1eacc782/?p=WXSMCKu6XGqj3_gUcV2aEw

**Задания:** 9–10

---

## 🗂️ Структура репозитория

├── Тема 1/

│ ├── ДВ_3_ИИ_Лекция_1.pptx

│ ├── ДВ_3_ИИ_Практика_1-1.pptx

│ ├── ДВ_3_ИИ_Практика_1-2.pptx

│ ├── ДВ_3_ИИ_Практика_1-3.pptx

│ ├── Практика_1-1.R

│ ├── Практика_1-2.R

│ └── Практика_1-3.R

├── Тема 2/

│ ├── ДВ_3_ИИ_Лекция_2.pptx

│ ├── ДВ_3_ИИ_Практика_2-1.pptx

│ ├── ДВ_3_ИИ_Практика_2-2.pptx

│ ├── ДВ_3_ИИ_Практика_2-3.pptx

│ ├── Практика 2-1.R

│ ├── Практика 2-2.R

│ └── Практика 2-3.R

├── Тема 3/ 

│ ├── ДВ_3_ИИ_Лекция_3.pptx

│ ├── ДВ_3_ИИ_Практика_3-1.pptx

│ ├── ДВ_3_ИИ_Практика_3-2.pptx

│ ├── ДВ_3_ИИ_Практика_3-3.pptx

│ ├── Практика 3-1.r

│ ├── Практика 3-2.r

│ └── Практика 3-3.r

├── Тема 4/ # в разработке

├── Тема 5/

│ └── ДВ_3_ИИ_Лекция_5.pptx

└── DV_3_Zadaniia.docx # все задания курса


---

Репозиторий пополняется по мере прохождения курса.

Email для связи: r.mukhamadeeva@yandex.ru
