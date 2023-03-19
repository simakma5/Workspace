/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2023 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */

/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */
#define setbit(reg, bit)      ((reg) |= (1U << (bit)))
#define clearbit(reg, bit)    ((reg) &= (~(1U << (bit))))
#define togglebit(reg, bit)   ((reg) ^= (1U << (bit)))
#define getbit(reg, bit)      (((reg) & (1U << (bit))) >> (bit))
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
UART_HandleTypeDef huart2;

/* USER CODE BEGIN PV */
const uint16_t BCDto7seg[16] = {
	0b0000000000000011,	// 0
	0b0000000010011111,	// 1
	0b0000000000100101,	// 2
	0b0000000000001101,	// 3
	0b0000000010011001,	// 4
	0b0000000001001001,	// 5
	0b0000000001000001,	// 6
	0b0000000000011111,	// 7
	0b0000000000000001,	// 8
	0b0000000000001001,	// 9
	0b1111111111111111,	// the rest is to avoid undefined behaviour
	0b1111111111111111,
	0b1111111111111111,
	0b1111111111111111,
	0b1111111111111111,
	0b1111111111111111
};
/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_USART2_UART_Init(void);
/* USER CODE BEGIN PFP */
void EXTI0_IRQHandler();
void EXTI4_IRQHandler();
void TIM3_IRQHandler();
/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{
  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_USART2_UART_Init();
  /* USER CODE BEGIN 2 */

	/* display */
	// PA9: data
	// PA8: clock
	// PB5: strobe

  	/* buttons */
	// PA1: lap
	// PA4: start-stop
	// PB0: reset

	/* timers */
    // TIM2: stopwatch value counter
	// TIM3: 7seg display clock
	// TIM5: start-stop button debouncing clock

  	/* allow clock for */
	// GPIOA and GPIOB
	RCC->AHB1ENR |= (0b1 << 0) | (0b1 << 1);
	// TIM2, TIM3 and TIM5
	RCC->APB1ENR |= (0b1 << 0) | (0b1 << 1) | (0b1 << 3);
	// SYSCFG
	RCC->APB2ENR |= (0b1 << 14);

	/* GPIOA init */
	// set PA8, PA9, PA5 as digital output, PA4 as digital input, PA1 as an alternate function (TIM2_CH2)
	GPIOA->MODER &= ~((0b11 << 2) | (0b11 << 8) | (0b11 << 10) | (0b11 << 16) | (0b11 << 18));
	GPIOA->MODER |= (0b10 << 2) | (0b01 << 10) | (0b01 << 16) | (0b01 << 18);
	// set low speed for all pins
	GPIOA->OSPEEDR &= ~((0b11 << 2) | (0b11 << 8) | (0b11 << 10) | (0b11 << 16) | (0b11 << 18));
	// disable PU/PD for all pins
	GPIOA->PUPDR &= ~((0b11 << 2) | (0b11 << 8) | (0b11 << 10) | (0b11 << 16) | (0b11 << 18));
	// PA1 alternate function = TIM2_CH2
	GPIOA->AFR[0] |= (0b0001 << 4);	// PA1 AF = TIM2_CH2

	/* EXTI4 init (PA4) */
	// allow interrupts from EXTI4
	EXTI->IMR |= (0b1 << 4);
	// set falling-edge trigger
	EXTI->FTSR |= (0b1 << 4);
	// clear the choice of GPIOx for EXTI4
	SYSCFG->EXTICR[1] &= ~(0b1111 << 0);
	// choose GPIOA for EXTI4
	SYSCFG->EXTICR[1] |= (0b0000 << 0);
	// enable interrupt EXTI4 in the NVIC controller
	NVIC_EnableIRQ(EXTI4_IRQn);

	/* GPIOB init */
	// set PB0 as a digital input and PB5 as a digital output
	GPIOB->MODER &= ~((0b11 << 0) | (0b11 << 10));
	GPIOB->MODER |= (0b01 << 10);
	// set low speed for both pins
	GPIOB->OSPEEDR &= ~((0b11 << 0) | (0b11 << 10));
	// disable PU/PD for both pins
	GPIOB->PUPDR &= ~((0b11 << 0) | (0b11 << 10));

	/* EXTI0 init (PB0) */
	// allow interrupts from EXTI0
	EXTI->IMR |= (0b1 << 0);
	// set falling-edge trigger
	EXTI->FTSR |= (0b1 << 0);
	// clear the choice of GPIOx for EXTI0
	SYSCFG->EXTICR[0] &= ~(0b1111 << 0);
	// choose GPIOB for EXTI0
	SYSCFG->EXTICR[0] |= (0b0001 << 0);
	// enable interrupt EXTI0 in the NVIC controller
	NVIC_EnableIRQ(EXTI0_IRQn);

	/* TIM2 init (stopwatch value) */
	// set prescaler value f_TIM2 = 1600 Hz = 2^4 * 100 Hz, where 2^4 is the value of the counter, i.e. T_TIM2 = 10 ms
	TIM2->PSC = (52500-1);
	TIM2->CNT = 0;
	// update counter registers
	TIM2->EGR |= (0b1 << 0);

	// set CC2 as input, IC2 as TI2
	TIM2->CCMR1 |= (0b01 << 8);
	// enable filtering, f_sample = f_TIM2, N = 8
	TIM2->CCMR1 |= (0b0011 << 12);
	// capture when LOW/FALLING EDGE
	TIM2->CCER |= (0b1 << 5);
	// enable capture/compare on channel 2
	TIM2->CCER |= (0b1 << 4);

	/* TIM3 init (strobe) */
	// enable update interrupt
	TIM3->DIER |= (0b1 << 0);
	NVIC_EnableIRQ(TIM3_IRQn);
	// TIM3 prescaler (f = 84 MHz)
	TIM3->PSC = 0;
	// automatic reload with f_TIM3 = (84*10^6)/([PSC+1]*[ARR+1]) Hz = 168kHz for PSC = 0 and ARR = 499
	TIM3->ARR = 499;
	// enable TIM3
	TIM3->CR1 |= (0b1 << 0);

	/* TIM5 init (debouncer) */
	// f_TIM5 = 1 kHz, i.e. T_TIM5 = 1 ms
	TIM5->PSC = (84000-1);
	TIM5->CNT = 0;
	// update counter registers
	TIM5->EGR |= (0b1 << 0);
	// enable TIM5
	TIM5->CR1 |= (0b1 << 0);

  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  while (1)
  {
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure the main internal regulator output voltage
  */
  __HAL_RCC_PWR_CLK_ENABLE();
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE2);

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState = RCC_HSI_ON;
  RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSI;
  RCC_OscInitStruct.PLL.PLLM = 16;
  RCC_OscInitStruct.PLL.PLLN = 336;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV4;
  RCC_OscInitStruct.PLL.PLLQ = 7;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_2) != HAL_OK)
  {
    Error_Handler();
  }
}

/**
  * @brief USART2 Initialization Function
  * @param None
  * @retval None
  */
static void MX_USART2_UART_Init(void)
{

  /* USER CODE BEGIN USART2_Init 0 */

  /* USER CODE END USART2_Init 0 */

  /* USER CODE BEGIN USART2_Init 1 */

  /* USER CODE END USART2_Init 1 */
  huart2.Instance = USART2;
  huart2.Init.BaudRate = 115200;
  huart2.Init.WordLength = UART_WORDLENGTH_8B;
  huart2.Init.StopBits = UART_STOPBITS_1;
  huart2.Init.Parity = UART_PARITY_NONE;
  huart2.Init.Mode = UART_MODE_TX_RX;
  huart2.Init.HwFlowCtl = UART_HWCONTROL_NONE;
  huart2.Init.OverSampling = UART_OVERSAMPLING_16;
  if (HAL_UART_Init(&huart2) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN USART2_Init 2 */

  /* USER CODE END USART2_Init 2 */

}

/**
  * @brief GPIO Initialization Function
  * @param None
  * @retval None
  */
static void MX_GPIO_Init(void)
{
  GPIO_InitTypeDef GPIO_InitStruct = {0};

  /* GPIO Ports Clock Enable */
  __HAL_RCC_GPIOC_CLK_ENABLE();
  __HAL_RCC_GPIOH_CLK_ENABLE();
  __HAL_RCC_GPIOA_CLK_ENABLE();
  __HAL_RCC_GPIOB_CLK_ENABLE();

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin : B1_Pin */
  GPIO_InitStruct.Pin = B1_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_FALLING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(B1_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : LD2_Pin */
  GPIO_InitStruct.Pin = LD2_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(LD2_GPIO_Port, &GPIO_InitStruct);

}

/* USER CODE BEGIN 4 */
void EXTI0_IRQHandler(){
	// disable TIM2 counter
	TIM2->CR1 &= ~(0b1 << 0);
	// clear TIM2 value
	TIM2->CNT = 0;
	// disable capture/compare on channel 2
	TIM2->CCER &= ~(0b1 << 4);
	// set CC2 as output so it can be reset
	TIM2->CCMR1 &= ~(0b11 << 8);
	// clear capture/compare on channel 2
	TIM2->CCR2 = 0;
	// set CC2 as input, IC2 as TI2
	TIM2->CCMR1 |= (0b01 << 8);
	// update counter registers
	TIM2->EGR |= (0b1 << 0);
	// enable capture/compare on channel 2
	TIM2->CCER |= (0b1 << 4);
	// clear EXTI0 interrupt flag
	EXTI->PR |= (0b1 << 0);
}

void EXTI4_IRQHandler(){
	// clear EXTI4 interrupt flag
	EXTI->PR |= (0b1 << 4);

	/* debouncing */
	// clear previous value
	volatile static int previous = 0;
	// load current value of TIM5
	volatile int current = TIM5->CNT;
	// prevent repeated activation
	if((current - previous) < 700) return;
	// prepare debouncer for next call
	previous = current;

	// toggle TIM2 counter enable (displayed count)
	togglebit(TIM2->CR1, 0);
}

void TIM3_IRQHandler(){
	static int digit = 0;
	static int count = 0;
	static uint8_t cislice[4] = {0};
	static uint16_t data = 0;
	if(count == 0) {
		if(digit == 0) {
			// check whether capture/compare occured
			uint32_t tmp;
			if(TIM2->CCR2 == 0) {
				// TIM2->CNT divided by 2^4, i.e. last digit = 0.01 s
				tmp = TIM2->CNT >> 4;
			} else {
				// TIM2->CC2 divided by 2^4, i.e. last digit = 0.01 s
				tmp = TIM2->CCR2 >> 4;
			}

			for(int i=3; i>=0; i--) {
				cislice[i] = tmp % 10;
				tmp = tmp / 10;
			}
		}
		// data = 1234xxxxABCDEFGH
		if(digit == 0 && cislice[0] == 0) {
			// disable tens of seconds when unused
			data = 0x00FF;
		} else if(digit != 1) {
			data = (0b1 << (15-digit)) | BCDto7seg[cislice[digit]];
		} else {
			// enable H segment
			data = (0b1 << (15-digit)) | (BCDto7seg[cislice[digit]] & 0b11111110);
		}
	}
	// write bit to DATA for 7seg
	uint16_t bit = (data >> count) & 0b1; // & 0b1 just for safety
	if(bit != 0) {
		GPIOA->BSRR |= (0b1 << 9);	// PA9='1'
	} else {
		GPIOA->BSRR |= (0b1 << 25);	// PA9='0'
	}
	GPIOA->BSRR |= (0b1 << 8);	// clock ↑
	__NOP();
	__NOP();
	__NOP();
	__NOP();
	__NOP();
	GPIOA->BSRR |= (0b1 << 24);	// clock ↓

	if(count == 15) {
		count = 0;
		// STROBE (PB5)
		GPIOB->BSRR |= (0b1 << 5);
		__NOP();
		__NOP();
		__NOP();
		__NOP();
		__NOP();
		GPIOB->BSRR |= (0b1 << 21);

		if(digit == 3) {
			digit = 0;
		} else {
			digit++;
		}
	} else {
		count++;
	}

	// clear update interrupt flag
	TIM3->SR &= ~(0b1 << 0);
}
/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
