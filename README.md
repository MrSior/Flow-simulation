# Flow-simulation

## Реализованно
- Ввод из файла
- Шаблонный ```Fixed``` в файле ```Fixed.h```
- Шаблонный класс симулятор в файле ```FluidEngine.h```
- Реализовано преобразование данных из рантайма в статический набор параметров для класса симулятора
- Возможность сохранить параметры симуляции с возможность запуска с данной точки

## Сборка проекта
```shell
mkdir cmake-build
cd cmake-build
cmake ..
cmake --build .
```
## Запуск и использованиk
### Пример:
```shell
./FluidEngine --p-type="FIXED(32,9)" --v-type="FIXED(32,9)" --v-flow-type="FIXED(32,9)" --input-file=../input.txt --save-file=../SimSave.txt
```
### Опции
- Допустимые типы для переменных ```--p-type```, ```--v-type```, ```--v-flow-type``` указаны в ```CmakeLists.txt``` в строчке ```TYPES``` внтури блока
    ```target_compile_definitions()```. Там же в ```SIZES``` указываются какие статические размеры поля надо скомпелировать.

- В параметре ```--input-file``` указывается файл с стартовой конфигурацией
- В параметре ```--save-file``` указывается файл для сохранения
- **Сохранение** осуществляется путём отправки сигнала ```SIGQUIT``` (сочетание клавиш в терминале ```CTR + 4```)
