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
#define setbit(reg, bit)		((reg) |= (1U << (bit)))
#define clearbit(reg, bit)		((reg) &= (~(1U << (bit))))
#define togglebit(reg, bit)		((reg) ^= (1U << (bit)))
#define getbit(reg, bit) 		(((reg) & (1U << (bit))) >> (bit))

#define P1 15000-1
#define P2 5000-1
#define P3 1000-1
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
UART_HandleTypeDef huart2;

/* USER CODE BEGIN PV */

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_USART2_UART_Init(void);
/* USER CODE BEGIN PFP */
void ledConfig();
void timerConfig();
void buttonConfig();
void TIM2_IRQHandler();
void EXTI15_10_IRQHandler();
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
  // MX_USART2_UART_Init();
  /* USER CODE BEGIN 2 */
  ledConfig();
  timerConfig();
  buttonConfig();
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
void ledConfig(){
  /* LED PA5 init */
  // allow clock for GPIOA
  setbit(RCC->AHB1ENR, 0);
  // reset GPIOA->MODER for pin 5
  clearbit(GPIOA->MODER, 11);
  clearbit(GPIOA->MODER, 10);
  // set GPIOA->MODER for pin 5 to 01 (output)
  setbit(GPIOA->MODER, 10);
  // set push-pull
  clearbit(GPIOA->OTYPER, 5);
  // disable pull-up/pull-down (GPIOA->PUPDR &= ~(0b11 << 10);)
  clearbit(GPIOA->PUPDR, 11);
  clearbit(GPIOA->PUPDR, 10);
}

void timerConfig(){
  /* TIM2 init */
  // allow clock
  setbit(RCC->APB1ENR, 0);
  // allow TIM2 interrupt
  setbit(TIM2->DIER, 0);
  NVIC_EnableIRQ(TIM2_IRQn);
  // prescaler setting
  TIM2->PSC = 4199;		// sets scaling ratio to (4199+1),
  	  	  	  	  	  	// effective clock clk/4200,
  	  	  	  	  	  	// results in a 50 us step
  // default auto-reload value
  TIM2->ARR = P1;
  // default timer register value
  TIM2->CNT = P1;
  // set to count down to zero
  setbit(TIM2->CR1, 4);
  // turn counter on
  setbit(TIM2->CR1, 0);
}

void buttonConfig(){
  /* Button PC13 init */
  // allow clock for GPIOC
  setbit(RCC->AHB1ENR, 2);
  // set PC13 as a digital input (GPIOC->MODER &= ~(0b11 << 26);)
  clearbit(GPIOC->MODER, 27);
  clearbit(GPIOC->MODER, 26);
  // disable PU/PD on PC13 (GPIOC->PUPDR &= ~(0b11 << 26);)
  clearbit(GPIOC->PUPDR, 27);
  clearbit(GPIOC->PUPDR, 26);

  /* EXTI13 init */
  // allow clock for SYSCFG (system configuration controller)
  setbit(RCC->APB2ENR, 14);
  // allow interrupts from EXTI13
  setbit(EXTI->IMR, 13);
  // select falling-edge trigger
  setbit(EXTI->FTSR, 13);
  // clear EXTI13 configuration
  SYSCFG->EXTICR[3] &= ~(0b1111 << 4);
  // set GPIOC as trigger for EXTI13
  SYSCFG->EXTICR[3] |= (0b0010 << 4);
  // allow EXTI15-10 interrupts (1 vector)
  NVIC_EnableIRQ(EXTI15_10_IRQn);
}

void TIM2_IRQHandler(){
  // clear update interrupt flag
  clearbit(TIM2->SR, 0);
  // toggle the state of LED PA5
  togglebit(GPIOA->ODR, 5);
}

void EXTI15_10_IRQHandler(){
  // clear EXTI13 interrupt flag
  setbit(EXTI->PR, 13);
  // load the value ARR of the counter TIM2
  uint32_t temp = TIM2->ARR;
  // set the next value of TIM2->ARR
  if(temp == P1){
	  TIM2->ARR = P2;
  } else if(temp == P2) {
  	TIM2->ARR = P3;
  } else {
  	TIM2->ARR = P1;
  }
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
