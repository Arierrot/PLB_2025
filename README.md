# Portuguese League of Bioinformatics 2025 - Gene Miners

Este repositorio contiene las soluciones a los cinco desafíos de la **Portuguese League of Bioinformatics 2025**, desarrolladas en pareja bajo el equipo **Gene Miners**. 

## Descripción

Los desafíos cubren diferentes áreas de la bioinformática, incluyendo estimación de tasas de mutación, recuperación de datos genómicos, análisis de anotaciones GO, identificación de mutaciones patogénicas y análisis de aptámeros.

## Desafíos y Soluciones

### 1️⃣ **Estimación de tasa de mutación con bootstrapping** (`challenge_I.py`)
Calcula la tasa de mutación por sitio y generación a partir de un archivo de entrada, aplicando bootstrapping para obtener intervalos de confianza.

### 2️⃣ **Extracción de genes desde NCBI Entrez API** (`challenge_II.py`)
Consulta la base de datos de NCBI para recuperar una lista de genes asociados a un organismo dado.

### 3️⃣ **Análisis de anotaciones GO de UniProt** (`challenge_III.py`)
Procesa un archivo de anotación GO para extraer y clasificar los 10 términos GO más relevantes para cada categoría (Proceso Biológico, Función Molecular, Componente Celular).

### 4️⃣ **Identificación de mutaciones exclusivas en cepas patógenas** (`challenge_IV.py`)
Compara secuencias de proteínas alineadas entre cepas patógenas y no patógenas para detectar mutaciones exclusivas.

### 5️⃣ **Análisis de aptámeros** (`challenge_V.py`)
Calcula el contenido GC, la longitud y la subcadena común más larga entre dos secuencias de aptámeros.

## Notas
- Los archivos de entrada deben seguir el formato indicado en cada desafío.
- Se recomienda ejecutar los scripts desde el directorio donde se encuentran los archivos de entrada.
- La extracción de datos desde NCBI puede tomar tiempo dependiendo del tamaño del organismo consultado.

---

**Equipo: Gene Miners 🚀**
Antía Dorado Valín | Pablo Torreira Pardo
