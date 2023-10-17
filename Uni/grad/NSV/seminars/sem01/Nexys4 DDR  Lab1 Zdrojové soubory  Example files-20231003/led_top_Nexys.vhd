----------------------------------------------------------------------------------
-- Company: CTU Prague
-- Engineer: Pavel Hazdra
-- 
-- Create Date: 12.07.2018 16:50:38
-- Design Name: 
-- Module Name: led_top - Behavioral
-- Project Name: 
-- Target Devices: 
-- Tool Versions: 
-- Description: 
-- 
-- Dependencies: 
-- 
-- Revision:
-- Revision 0.01 - File Created
-- Additional Comments:
-- 
----------------------------------------------------------------------------------


library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.STD_LOGIC_ARITH.ALL;
use IEEE.STD_LOGIC_UNSIGNED.ALL;


-- Uncomment the following library declaration if using
-- arithmetic functions with Signed or Unsigned values
--use IEEE.NUMERIC_STD.ALL;

-- Uncomment the following library declaration if instantiating
-- any Xilinx leaf cells in this code.
--library UNISIM;
--use UNISIM.VComponents.all;

entity led_top is
    Port ( clk_i : in STD_LOGIC;
           rstn_i : in STD_LOGIC;
           led_o : out STD_LOGIC_VECTOR (15 downto 0);
           disp_seg_o : out STD_LOGIC_VECTOR (7 downto 0);
           disp_an_o : out STD_LOGIC_VECTOR (7 downto 0));
end led_top;

architecture Behavioral of led_top is


component seven_seg is
   port (
        hex :  in  std_logic_vector(3 downto 0);
       ss  :  out  std_logic_vector(7 downto 0)
     );
end component seven_seg;


signal cnt : std_logic_vector (31 downto 0);

begin

seven_seg_inst: seven_seg
   port map (
       hex(3 downto 0) => cnt(27 downto 24),
       ss(7 downto 0)  => disp_seg_o 
    );

process (clk_i, rstn_i)
begin

    if rstn_i = '0' then cnt <= (others => '0');
        elsif clk_i = '1' and clk_i'event then
            cnt <= cnt + 1;
    end if;

end process;

disp_an_o <= "00001111";
led_o <= cnt(31 downto 16);

end Behavioral;
